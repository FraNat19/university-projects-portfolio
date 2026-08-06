#!/usr/bin/env python3

"""

LAYER 3: INDEXES + VERIFY



- Non crea constraints (le crea già layer3_fixed.py con ensure_schema()).

- Crea solo indici utili.

- Provide verify per contare nodi layer3.

"""



import os

import argparse

from neo4j import GraphDatabase



NEO4J_URI = os.getenv("NEO4J_URI", "bolt://127.0.0.1:7687")

NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")

NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "")



driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))



def create_layer3_indexes():

    print("\n📊 LAYER 3: creating indexes")

    print("=" * 70)



    indexes = [

        "CREATE INDEX RiskCategory_severity IF NOT EXISTS FOR (r:RiskCategory) ON (r.severity_level)",

        "CREATE INDEX RiskCategory_frequency IF NOT EXISTS FOR (r:RiskCategory) ON (r.frequency)",

        "CREATE INDEX ProfessionalTag_sector IF NOT EXISTS FOR (p:ProfessionalTag) ON (p.sector)",

        "CREATE INDEX ProfessionalTag_frequency IF NOT EXISTS FOR (p:ProfessionalTag) ON (p.frequency)",

        "CREATE INDEX Topic_frequency IF NOT EXISTS FOR (t:Topic) ON (t.frequency)",

    ]



    with driver.session() as session:

        for stmt in indexes:

            session.run(stmt)

            print(f" ✅ {stmt.split()[2]}")



    print("\n✅ Indexes created.\n")



def verify_layer3():

    print("\n🔍 LAYER 3: verify")

    print("=" * 70)



    with driver.session() as session:

        rec = session.run("""

        RETURN

          count{(r:RiskCategory)} AS risks,

          count{(p:ProfessionalTag)} AS professions,

          count{(t:Topic)} AS topics,

          count{(:TechnicalDocument)-[:ADDRESSES]->(:RiskCategory)} AS doc_risk_links,

          count{(:TechnicalDocument)-[:TARGETS]->(:ProfessionalTag)} AS doc_prof_links,

          count{(:TechnicalDocument)-[:DISCUSSES]->(:Topic)} AS doc_topic_links,

          count{(:RiskCategory)-[:GOVERNED_BY]->(:CanonicalLaw)} AS risk_law_links,

          count{(:ProfessionalTag)-[:REGULATED_BY]->(:CanonicalLaw)} AS prof_law_links

        """).single()



    print(f"RiskCategory nodes:      {rec['risks']:,}")

    print(f"ProfessionalTag nodes:   {rec['professions']:,}")

    print(f"Topic nodes:             {rec['topics']:,}")

    print(f"Doc-ADDRESSES-Risk:      {rec['doc_risk_links']:,}")

    print(f"Doc-TARGETS-Profession:  {rec['doc_prof_links']:,}")

    print(f"Doc-DISCUSSES-Topic:     {rec['doc_topic_links']:,}")

    print(f"Risk-GOVERNED_BY-Law:    {rec['risk_law_links']:,}")

    print(f"Prof-REGULATED_BY-Law:   {rec['prof_law_links']:,}")

    print("\n✅ Verify done.\n")



def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--create", action="store_true", help="Create indexes")

    parser.add_argument("--verify", action="store_true", help="Verify layer3 content")

    args = parser.parse_args()



    if not args.create and not args.verify:

        create_layer3_indexes()

        verify_layer3()

        return



    if args.create:

        create_layer3_indexes()

    if args.verify:

        verify_layer3()



if __name__ == "__main__":

    try:

        main()

    finally:

        driver.close()

