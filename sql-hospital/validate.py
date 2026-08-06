"""Run schema.sql, seed_data.sql and queries.sql without a MySQL server.

Translates the MySQL dialect to SQLite in memory, loads the schema and the
sample records, then executes every validation query and prints its result.
Foreign keys are enforced, so a broken reference fails the run.

Francesco Natali
"""

import re
import sqlite3
import sys
from datetime import date
from pathlib import Path

HERE = Path(__file__).parent


def years_between(start, end):
    if start is None or end is None:
        return None
    a = date.fromisoformat(str(start)[:10])
    b = date.fromisoformat(str(end)[:10])
    return b.year - a.year - ((b.month, b.day) < (a.month, a.day))


def to_sqlite(sql):
    sql = re.sub(r"^\s*(DROP DATABASE|CREATE DATABASE|USE)\b.*?;", "", sql,
                 flags=re.I | re.M | re.S)
    sql = re.sub(r"\)\s*ENGINE=InnoDB", ")", sql, flags=re.I)
    sql = re.sub(r"\bINT\s+AUTO_INCREMENT\s+PRIMARY KEY\b", "INTEGER PRIMARY KEY AUTOINCREMENT",
                 sql, flags=re.I)
    sql = re.sub(r"\bENUM\s*\(([^)]*)\)", lambda m: f"TEXT CHECK (1 OR value IN ({m.group(1)}))",
                 sql, flags=re.I)
    sql = re.sub(r"TEXT CHECK \(1 OR value IN \(([^)]*)\)\)", "TEXT", sql)
    sql = re.sub(r"\b(TINYINT|SMALLINT|MEDIUMINT|BIGINT)\b", "INTEGER", sql, flags=re.I)
    sql = re.sub(r"\bBOOLEAN\b", "INTEGER", sql, flags=re.I)
    sql = re.sub(r"\bDATETIME\b|\bDATE\b|\bTIME\b", "TEXT", sql, flags=re.I)
    sql = re.sub(r"\bCHAR\s*\(\d+\)|\bVARCHAR\s*\(\d+\)", "TEXT", sql, flags=re.I)
    sql = re.sub(r"\bTRUE\b", "1", sql, flags=re.I)
    sql = re.sub(r"\bFALSE\b", "0", sql, flags=re.I)
    sql = re.sub(r"TIMESTAMPDIFF\s*\(\s*YEAR\s*,", "YEARS_BETWEEN(", sql, flags=re.I)
    return sql


def statements(sql):
    for raw in sql.split(";"):
        s = "\n".join(l for l in raw.splitlines() if not l.strip().startswith("--")).strip()
        if s:
            yield s


def main():
    con = sqlite3.connect(":memory:")
    con.create_function("YEARS_BETWEEN", 2, years_between)
    con.execute("PRAGMA foreign_keys = ON")

    for name in ("schema.sql", "seed_data.sql"):
        text = to_sqlite((HERE / name).read_text(encoding="utf-8"))
        count = 0
        for stmt in statements(text):
            try:
                con.execute(stmt)
                count += 1
            except sqlite3.Error as exc:
                print(f"FAIL in {name}: {exc}\n{stmt[:200]}")
                return 1
        print(f"{name}: {count} statements executed")

    tables = [r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name")]
    print(f"\ntables created: {len(tables)}")
    for t in tables:
        n = con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        print(f"  {t:22} {n:4} rows")

    broken = con.execute("PRAGMA foreign_key_check").fetchall()
    print(f"\nforeign key violations: {len(broken)}")
    if broken:
        for row in broken[:10]:
            print(" ", row)
        return 1

    text = to_sqlite((HERE / "queries.sql").read_text(encoding="utf-8"))
    blocks = re.split(r"\n(?=-- Q\d)", (HERE / "queries.sql").read_text(encoding="utf-8"))
    print("\n" + "=" * 66)
    for block in blocks:
        m = re.match(r"-- (Q\d) - (.*)", block.strip())
        if not m:
            continue
        label, title = m.groups()
        body = to_sqlite(block)
        for stmt in statements(body):
            cur = con.execute(stmt)
            if stmt.lstrip().upper().startswith("SELECT"):
                rows = cur.fetchall()
                cols = [d[0] for d in cur.description]
                print(f"\n{label} - {title}")
                print("  " + " | ".join(cols))
                for r in rows[:6]:
                    print("  " + " | ".join("" if v is None else str(v) for v in r))
                if len(rows) > 6:
                    print(f"  ... {len(rows) - 6} more rows")
                if not rows:
                    print("  (no rows)")
            else:
                con.commit()
                print(f"\n{label} - {title}\n  {cur.rowcount} row(s) affected")

    print("\n" + "=" * 66)
    print("all statements executed, referential integrity holds")
    return 0


if __name__ == "__main__":
    sys.exit(main())
