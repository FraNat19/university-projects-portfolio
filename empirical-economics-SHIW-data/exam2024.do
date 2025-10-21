clear

cd "X:\APPLIEDECONOMICS\EMPEC24\SHIW10\"

use "q10c2.dta"


gen SCONTO = .

* Half of the sample asked with 10% first
replace SCONTO = 10 if sconto1 == 1
replace SCONTO = 20 if sconto1 == 1 & sconto2 == 1
replace SCONTO = 5 if sconto1 == 2 & sconto3 == 1
replace SCONTO = 2 if sconto1 == 2 & sconto2 == 2 & sconto3 == 2 & sconto4 == 1

* Half of the sample asked with 20% first
replace SCONTO = 20 if sconto21 == 1
replace SCONTO = 10 if sconto21 == 2 & sconto22 == 1
replace SCONTO = 5 if sconto21 == 2 & sconto22 == 2 & sconto23 == 1
replace SCONTO = 2 if sconto21 == 2 & sconto22 == 2 & sconto23 == 2 & sconto24 == 1

replace SCONTO = 0 if sconto1 == 2 & sconto2 == 2 & sconto3 == 2 & sconto4 == 2
replace SCONTO = 0 if sconto21 == 2 & sconto22 == 2 & sconto23 == 2 & sconto24 == 2

tab SCONTO



graph bar, over(SCONTO)



joinby nquest using "carcom10.dta"

joinby nquest using "q10a.dta"

joinby nquest using "lavoro.dta"

joinby nquest using "allb1.dta"

joinby nquest using "q10d.dta"

joinby nquest using "q10e.dta"

keep if nord==1

svyset nquest [pweight=pesofit]

tab SCONTO


recode sex (1=0) (2=1)

svy: tab SCONTO sex, row

tab SCONTO sex, column
tab SCONTO sex, row

svy: oprobit SCONTO sex
svy: oprobit SCONTO sex i.studio


gen hstud=0
replace hstud=1 if studio<5
replace hstud=2 if studio==5
replace hstud=3 if studio>5

svy: oprobit SCONTO sex i.hstud

svy: oprobit SCONTO sex i.hstud i.stupcf i.stumcf

gen fed=(stupcf>=3&stupcf<=6)

gen med=(stumcf>=3&stumcf<=6)

svy: oprobit SCONTO sex i.hstud fed med


svy: tab SCONTO if eta<30
svy: tab SCONTO if eta>=30&eta<40
svy: tab SCONTO if eta>=40&eta<50
svy: tab SCONTO if eta>=50&eta<60
svy: tab SCONTO if eta>=60&eta<70
svy: tab SCONTO if eta>=70&eta<80
svy: tab SCONTO if eta>=80


svy: oprobit SCONTO sex i.hstud eta 

gen risk=0
replace risk=1 if risfin==1
replace risk=2 if risfin==2
replace risk=3 if risfin==3
replace risk=4 if risfin==4



svy: oprobit SCONTO sex i.hstud eta risk

gen income=0
replace income=1 if ylm<=10000
replace income=2 if ylm>10000 & ylm<=20000
replace income=3 if ylm>20000 & ylm<=30000
replace income=4 if ylm>30000



svy: oprobit SCONTO sex i.hstud eta risk i.income


gen rispneg = -rispbass

replace rispneg = 0 if rispneg == .

replace rispalt = 0 if rispalt == .

gen RISPARMI = rispalt + rispneg

tab RISPARMI

gen Risparmi=0
replace Risparmi=1 if RISPARMI>= -50000
replace Risparmi=2 if RISPARMI> -50000 & Risparmi<= -25000
replace Risparmi=3 if RISPARMI> -25000 & Risparmi<0
replace Risparmi=4 if RISPARMI==0
replace Risparmi=5 if RISPARMI>0 & Risparmi<25000
replace Risparmi=6 if RISPARMI>=25000 & Risparmi<50000
replace Risparmi=7 if RISPARMI>=50000


egen risp1 = cut(RISPARMI), at(-50000, -25000, -0.1, 0.1, 25000, 50000, 500000)

recode RISPARMI (-200000/-49999=1 ) (-50000/-24999=2 ) (-25000/-0.1=3 ) (0/0.1=4 ) (0.2/24999=5 ) (25000/49999=6 ) (50000/500000=5 ), gen(risp)


recode RISPARMI (-200000/-24999=1) (-25000/-0.1=2) (0/0.1=3) (0.2/24999=4) (25000/500000=5), gen(saving)




egen prec = cut(precauz), at(20000, 50000, 100000, 200000, 500000)

recode precauz (0/20000=1 ) (20001/50000=2 ) (50001/100000=3) (100001/200000=4 ) (200001/500000=5 ) (500001/2000000=6) , gen(precauzione)
recode precauz (0/20000=1 ) (20001/50000=2 ) (50001/100000=3) (100001/200000=4 ) (200001/2000000=5 ) , gen(precaution)


svy: oprobit SCONTO sex i.hstud eta risk i.income saving precaution


gen properties=0
replace properties=1 if altrab==2 & altrfab==2 & teragr==2 & ternagr==2
replace properties=2 if altrab==1 | altrfab==1
replace properties=3 if teragr==1 | ternagr==1
replace properties=4 if (altrab==1 | altrfab==1) & (teragr==1 | ternagr==1)


gen abit2=.
replace abit2=0 if altrab==2 & altrfab==2 & teragr==2 & ternagr==2
replace abit2=1 if altrab==1 | altrfab==1 | teragr==1 | ternagr==1

recode pianoris (2=0) (1=1)



replace varcons=0 if varcons==4 | varcons==5

replace condgen=1 if condgen==1 | condgen==2
replace condgen=2 if condgen==3
replace condgen=3 if condgen==4
replace condgen=4 if condgen==5 | condgen==6


svy: oprobit SCONTO sex i.hstud eta risk i.income saving precaution i.condgen



gen debt=0
replace debt=1 if debita1==2 & debita2==2
replace debt=2 if debita1==1 
replace debt=3 if debita2==1
replace debt=4 if debita1==1 & debita2==1
replace debt = 0 if debt == .

gen nfam=0
replace nfam=1 if ncomp==1 | ncomp==2
replace nfam=2 if ncomp==3
replace nfam=3 if ncomp==4
replace nfam=4 if ncomp==5
replace nfam=5 if ncomp>=6

rename varred var_income
recode var_income (5=4) 

rename q working_status

gen age_class=1 if eta<25
replace age_class=2 if eta>=25&eta<33
replace age_class=3 if eta>=33&eta<41
replace age_class=4 if eta>=41&eta<51
replace age_class=5 if eta>=51&eta<61
replace age_class=6 if eta>=61&eta<71
replace age_class=7 if eta>=71&eta<81
replace age_class=8 if eta>=81&eta<90
replace age_class=9 if eta>=90


oprobit SCONTO sex eta i.hstud i.area3 staciv i.risk i.income i.saving i.precaution i.condgen i.debt i.properties pianoris

gen id_group=.
replace id_group=0 if sconto1==1|sconto1==2|sconto2==1|sconto2==2|sconto3==1|sconto3==2|sconto4==1|sconto4==2
replace id_group=1 if sconto21==1|sconto21==2|sconto22==1|sconto22==2|sconto23==1|sconto23==2|sconto24==1|sconto24==2

oprobit SCONTO id_group sex eta i.hstud i.area3 staciv i.risk i.income saving i.precaution i.condgen i.debt i.properties pianoris 






