function dydt=SMC_mTOR_ODE(t,y,params) 
% SMC_mTOR_ODE.m 
 
% Assign names for parameters 
[rpar,tau,ymax,speciesNames]=params{:}; 
Akt = 1; 
AMPK = 2; 
AngII = 3; 
AngIIin = 4; 
Apoptosis = 5; 
AT1R = 6; 
AT2R = 7; 
betaCatenin = 8; 
calponin = 9; 
Col3a1 = 10; 
Collagen = 11; 
DegElastin = 12; 
Dsh = 13; 
Elastin = 14; 
Elk1 = 15; 
Eln = 16; 
Energy = 17; 
ERK12 = 18; 
FAK = 19; 
Fibrillin = 20; 
FOXO = 21; 
Frizzled = 22; 
Glucose = 23; 
GSK3 = 24; 
IGF = 25; 
IGFR = 26; 
IL6 = 27; 
IFNgamma = 28; 
Integrins = 29; 
IRS1 = 30; 
JAK = 31; 
JNK = 32; 
LAMP1 = 33; 
LAMP2 = 34; 
latTGFb = 35; 
Leucine = 36; 
MAPK = 37; 
MEK = 38; 
mitf = 39; 
MITF = 40; 
MitoMet = 41; 
MMP2 = 42; 
mTOR = 43; 
mTORC1 = 44; 
mTORC2 = 45; 
NFkB = 46; 
Oxygen = 47; 
p38 = 48; 
p4EBP1 = 49; 
PA = 50; 
PDGF = 51; 
PDGFR = 52; 
PDK1 = 53; 
PI3K = 54; 
PLD = 55; 
Raf = 56; 
Rag = 57; 
Ras = 58; 
RhoA = 59; 
RSK1 = 60; 
S6 = 61; 
S6K = 62; 
Shear = 63; 
SM22 = 64; 
SMA = 65; 
Smad23 = 66; 
SMMHC = 67; 
SOCs = 68; 
STAT = 69; 
Stress = 70; 
TGFbeta = 71; 
TGFR = 72; 
TIMP = 73; 
TSC12 = 74; 
Wnt = 75; 
Wnt5 = 76; 
dydt = zeros(76,1); 
dydt(Akt) = (OR(act(y(mTORC2),rpar(:,91)),act(y(PDK1),rpar(:,102)))*ymax(Akt) - y(Akt))/tau(Akt); 
dydt(AMPK) = (AND(rpar(:,19),act(y(AT1R),rpar(:,19)),inhib(y(Energy),rpar(:,19)),inhib(y(Oxygen),rpar(:,19)))*ymax(AMPK) - y(AMPK))/tau(AMPK); 
dydt(AngII) = (OR(act(y(AngIIin),rpar(:,41)),act(y(Stress),rpar(:,120)))*ymax(AngII) - y(AngII))/tau(AngII); 
dydt(AngIIin) = (rpar(1,1)*ymax(AngIIin) - y(AngIIin))/tau(AngIIin); 
dydt(Apoptosis) = (act(y(FOXO),rpar(:,66))*ymax(Apoptosis) - y(Apoptosis))/tau(Apoptosis); 
dydt(AT1R) = (act(y(AngII),rpar(:,39))*ymax(AT1R) - y(AT1R))/tau(AT1R); 
dydt(AT2R) = (act(y(AngII),rpar(:,40))*ymax(AT2R) - y(AT2R))/tau(AT2R); 
dydt(betaCatenin) = (inhib(y(GSK3),rpar(:,23))*ymax(betaCatenin) - y(betaCatenin))/tau(betaCatenin); 
dydt(calponin) = (OR(AND(rpar(:,15),inhib(y(Elk1),rpar(:,15)),act(y(RhoA),rpar(:,15))),OR(AND(rpar(:,30),act(y(Akt),rpar(:,30)),inhib(y(FOXO),rpar(:,30)),inhib(y(NFkB),rpar(:,30))),OR(AND(rpar(:,59),act(y(ERK12),rpar(:,59)),act(y(Shear),rpar(:,59))),AND(rpar(:,113),inhib(y(FOXO),rpar(:,113)),inhib(y(NFkB),rpar(:,113)),act(y(p38),rpar(:,113)),act(y(Smad23),rpar(:,113))))))*ymax(calponin) - y(calponin))/tau(calponin); 
dydt(Col3a1) = (OR(act(y(p38),rpar(:,95)),act(y(Smad23),rpar(:,116)))*ymax(Col3a1) - y(Col3a1))/tau(Col3a1); 
dydt(Collagen) = (OR(inhib(y(MMP2),rpar(:,24)),act(y(Col3a1),rpar(:,56)))*ymax(Collagen) - y(Collagen))/tau(Collagen); 
dydt(DegElastin) = (inhib(y(Elastin),rpar(:,14))*ymax(DegElastin) - y(DegElastin))/tau(DegElastin); 
dydt(Dsh) = (act(y(Frizzled),rpar(:,67))*ymax(Dsh) - y(Dsh))/tau(Dsh); 
dydt(Elastin) = (OR(inhib(y(MMP2),rpar(:,25)),AND(rpar(:,57),act(y(Eln),rpar(:,57)),act(y(Fibrillin),rpar(:,57))))*ymax(Elastin) - y(Elastin))/tau(Elastin); 
dydt(Elk1) = (act(y(ERK12),rpar(:,60))*ymax(Elk1) - y(Elk1))/tau(Elk1); 
dydt(Eln) = (OR(act(y(IGF),rpar(:,70)),act(y(Smad23),rpar(:,117)))*ymax(Eln) - y(Eln))/tau(Eln); 
dydt(Energy) = (rpar(1,2)*ymax(Energy) - y(Energy))/tau(Energy); 
dydt(ERK12) = (AND(rpar(:,83),inhib(y(AT2R),rpar(:,83)),act(y(MEK),rpar(:,83)))*ymax(ERK12) - y(ERK12))/tau(ERK12); 
dydt(FAK) = (act(y(Integrins),rpar(:,75))*ymax(FAK) - y(FAK))/tau(FAK); 
dydt(Fibrillin) = (rpar(1,3)*ymax(Fibrillin) - y(Fibrillin))/tau(Fibrillin); 
dydt(FOXO) = (inhib(y(Akt),rpar(:,11))*ymax(FOXO) - y(FOXO))/tau(FOXO); 
dydt(Frizzled) = (act(y(Wnt),rpar(:,131))*ymax(Frizzled) - y(Frizzled))/tau(Frizzled); 
dydt(Glucose) = (rpar(1,4)*ymax(Glucose) - y(Glucose))/tau(Glucose); 
dydt(GSK3) = (OR(AND(rpar(:,10),inhib(y(Akt),rpar(:,10)),inhib(y(Dsh),rpar(:,10))),AND(rpar(:,13),inhib(y(Dsh),rpar(:,13)),inhib(y(S6K),rpar(:,13))))*ymax(GSK3) - y(GSK3))/tau(GSK3); 
dydt(IGF) = (act(y(Stress),rpar(:,121))*ymax(IGF) - y(IGF))/tau(IGF); 
dydt(IGFR) = (act(y(IGF),rpar(:,71))*ymax(IGFR) - y(IGFR))/tau(IGFR); 
dydt(IL6) = (act(y(NFkB),rpar(:,93))*ymax(IL6) - y(IL6))/tau(IL6); 
dydt(IFNgamma) = (rpar(1,5)*ymax(IFNgamma) - y(IFNgamma))/tau(IFNgamma); 
dydt(Integrins) = (OR(act(y(Fibrillin),rpar(:,65)),act(y(Stress),rpar(:,122)))*ymax(Integrins) - y(Integrins))/tau(Integrins); 
dydt(IRS1) = (inhib(y(S6K),rpar(:,28))*ymax(IRS1) - y(IRS1))/tau(IRS1); 
dydt(JAK) = (OR(act(y(AT1R),rpar(:,43)),OR(AND(rpar(:,73),act(y(IL6),rpar(:,73)),inhib(y(SOCs),rpar(:,73))),act(y(IFNgamma),rpar(:,74))))*ymax(JAK) - y(JAK))/tau(JAK); 
dydt(JNK) = (OR(act(y(AMPK),rpar(:,37)),OR(act(y(AT1R),rpar(:,44)),OR(act(y(FAK),rpar(:,63)),OR(act(y(PDGFR),rpar(:,98)),act(y(TGFR),rpar(:,130))))))*ymax(JNK) - y(JNK))/tau(JNK); 
dydt(LAMP1) = (OR(act(y(betaCatenin),rpar(:,53)),act(y(MITF),rpar(:,85)))*ymax(LAMP1) - y(LAMP1))/tau(LAMP1); 
dydt(LAMP2) = (OR(act(y(betaCatenin),rpar(:,54)),act(y(MITF),rpar(:,86)))*ymax(LAMP2) - y(LAMP2))/tau(LAMP2); 
dydt(latTGFb) = (OR(act(y(Integrins),rpar(:,76)),act(y(MAPK),rpar(:,81)))*ymax(latTGFb) - y(latTGFb))/tau(latTGFb); 
dydt(Leucine) = (rpar(1,6)*ymax(Leucine) - y(Leucine))/tau(Leucine); 
dydt(MAPK) = (act(y(Raf),rpar(:,105))*ymax(MAPK) - y(MAPK))/tau(MAPK); 
dydt(MEK) = (act(y(MAPK),rpar(:,82))*ymax(MEK) - y(MEK))/tau(MEK); 
dydt(mitf) = (act(y(betaCatenin),rpar(:,55))*ymax(mitf) - y(mitf))/tau(mitf); 
dydt(MITF) = (AND(rpar(:,84),act(y(GSK3),rpar(:,84)),act(y(mitf),rpar(:,84)),act(y(mTORC1),rpar(:,84)))*ymax(MITF) - y(MITF))/tau(MITF); 
dydt(MitoMet) = (act(y(Glucose),rpar(:,68))*ymax(MitoMet) - y(MitoMet))/tau(MitoMet); 
dydt(MMP2) = (OR(AND(rpar(:,34),act(y(Akt),rpar(:,34)),act(y(JNK),rpar(:,34)),act(y(p38),rpar(:,34)),inhib(y(TIMP),rpar(:,34)),inhib(y(TSC12),rpar(:,34))),OR(AND(rpar(:,52),act(y(betaCatenin),rpar(:,52)),inhib(y(TIMP),rpar(:,52))),AND(rpar(:,58),act(y(ERK12),rpar(:,58)),inhib(y(TIMP),rpar(:,58)))))*ymax(MMP2) - y(MMP2))/tau(MMP2); 
dydt(mTOR) = (act(y(Akt),rpar(:,35))*ymax(mTOR) - y(mTOR))/tau(mTOR); 
dydt(mTORC1) = (OR(inhib(y(TSC12),rpar(:,29)),OR(AND(rpar(:,87),act(y(MitoMet),rpar(:,87)),act(y(mTOR),rpar(:,87))),OR(AND(rpar(:,96),act(y(mTOR),rpar(:,96)),act(y(PA),rpar(:,96))),AND(rpar(:,106),act(y(mTOR),rpar(:,106)),act(y(Rag),rpar(:,106))))))*ymax(mTORC1) - y(mTORC1))/tau(mTORC1); 
dydt(mTORC2) = (act(y(mTOR),rpar(:,88))*ymax(mTORC2) - y(mTORC2))/tau(mTORC2); 
dydt(NFkB) = (act(y(AT1R),rpar(:,45))*ymax(NFkB) - y(NFkB))/tau(NFkB); 
dydt(Oxygen) = (rpar(1,7)*ymax(Oxygen) - y(Oxygen))/tau(Oxygen); 
dydt(p38) = (OR(act(y(AT1R),rpar(:,46)),OR(act(y(PDGFR),rpar(:,99)),act(y(TGFR),rpar(:,127))))*ymax(p38) - y(p38))/tau(p38); 
dydt(p4EBP1) = (act(y(mTORC1),rpar(:,89))*ymax(p4EBP1) - y(p4EBP1))/tau(p4EBP1); 
dydt(PA) = (act(y(PLD),rpar(:,104))*ymax(PA) - y(PA))/tau(PA); 
dydt(PDGF) = (OR(act(y(Shear),rpar(:,110)),act(y(Stress),rpar(:,123)))*ymax(PDGF) - y(PDGF))/tau(PDGF); 
dydt(PDGFR) = (act(y(PDGF),rpar(:,97))*ymax(PDGFR) - y(PDGFR))/tau(PDGFR); 
dydt(PDK1) = (act(y(PI3K),rpar(:,103))*ymax(PDK1) - y(PDK1))/tau(PDK1); 
dydt(PI3K) = (OR(AND(rpar(:,42),act(y(AT1R),rpar(:,42)),act(y(JAK),rpar(:,42))),OR(AND(rpar(:,77),act(y(IGFR),rpar(:,77)),act(y(IRS1),rpar(:,77))),act(y(PDGFR),rpar(:,101))))*ymax(PI3K) - y(PI3K))/tau(PI3K); 
dydt(PLD) = (act(y(AT1R),rpar(:,47))*ymax(PLD) - y(PLD))/tau(PLD); 
dydt(Raf) = (act(y(Ras),rpar(:,107))*ymax(Raf) - y(Raf))/tau(Raf); 
dydt(Rag) = (act(y(Leucine),rpar(:,80))*ymax(Rag) - y(Rag))/tau(Rag); 
dydt(Ras) = (OR(act(y(AT1R),rpar(:,48)),OR(act(y(FAK),rpar(:,64)),OR(act(y(IGFR),rpar(:,72)),act(y(PDGFR),rpar(:,100)))))*ymax(Ras) - y(Ras))/tau(Ras); 
dydt(RhoA) = (OR(AND(rpar(:,92),inhib(y(AT2R),rpar(:,92)),act(y(mTORC2),rpar(:,92))),act(y(Ras),rpar(:,108)))*ymax(RhoA) - y(RhoA))/tau(RhoA); 
dydt(RSK1) = (act(y(ERK12),rpar(:,61))*ymax(RSK1) - y(RSK1))/tau(RSK1); 
dydt(S6) = (act(y(S6K),rpar(:,109))*ymax(S6) - y(S6))/tau(S6); 
dydt(S6K) = (act(y(mTORC1),rpar(:,90))*ymax(S6K) - y(S6K))/tau(S6K); 
dydt(Shear) = (rpar(1,8)*ymax(Shear) - y(Shear))/tau(Shear); 
dydt(SM22) = (OR(AND(rpar(:,16),inhib(y(Elk1),rpar(:,16)),act(y(RhoA),rpar(:,16))),OR(AND(rpar(:,20),inhib(y(ERK12),rpar(:,20)),act(y(Shear),rpar(:,20))),OR(AND(rpar(:,31),act(y(Akt),rpar(:,31)),inhib(y(FOXO),rpar(:,31)),inhib(y(NFkB),rpar(:,31))),AND(rpar(:,114),inhib(y(FOXO),rpar(:,114)),inhib(y(NFkB),rpar(:,114)),act(y(p38),rpar(:,114)),act(y(Smad23),rpar(:,114))))))*ymax(SM22) - y(SM22))/tau(SM22); 
dydt(SMA) = (OR(AND(rpar(:,17),inhib(y(Elk1),rpar(:,17)),act(y(RhoA),rpar(:,17))),OR(AND(rpar(:,21),inhib(y(ERK12),rpar(:,21)),act(y(Shear),rpar(:,21))),OR(AND(rpar(:,32),act(y(Akt),rpar(:,32)),inhib(y(FOXO),rpar(:,32)),inhib(y(NFkB),rpar(:,32))),AND(rpar(:,112),inhib(y(FOXO),rpar(:,112)),inhib(y(NFkB),rpar(:,112)),act(y(p38),rpar(:,112)),act(y(Smad23),rpar(:,112))))))*ymax(SMA) - y(SMA))/tau(SMA); 
dydt(Smad23) = (OR(act(y(AT1R),rpar(:,49)),AND(rpar(:,126),act(y(TGFR),rpar(:,126)),act(y(TSC12),rpar(:,126))))*ymax(Smad23) - y(Smad23))/tau(Smad23); 
dydt(SMMHC) = (OR(AND(rpar(:,18),inhib(y(Elk1),rpar(:,18)),act(y(RhoA),rpar(:,18))),OR(AND(rpar(:,33),act(y(Akt),rpar(:,33)),inhib(y(FOXO),rpar(:,33)),inhib(y(NFkB),rpar(:,33))),OR(act(y(ERK12),rpar(:,62)),AND(rpar(:,115),inhib(y(FOXO),rpar(:,115)),inhib(y(NFkB),rpar(:,115)),act(y(p38),rpar(:,115)),act(y(Smad23),rpar(:,115))))))*ymax(SMMHC) - y(SMMHC))/tau(SMMHC); 
dydt(SOCs) = (act(y(STAT),rpar(:,118))*ymax(SOCs) - y(SOCs))/tau(SOCs); 
dydt(STAT) = (act(y(JAK),rpar(:,78))*ymax(STAT) - y(STAT))/tau(STAT); 
dydt(Stress) = (rpar(1,9)*ymax(Stress) - y(Stress))/tau(Stress); 
dydt(TGFbeta) = (OR(AND(rpar(:,79),inhib(y(Fibrillin),rpar(:,79)),act(y(latTGFb),rpar(:,79))),OR(act(y(Shear),rpar(:,111)),act(y(Stress),rpar(:,124))))*ymax(TGFbeta) - y(TGFbeta))/tau(TGFbeta); 
dydt(TGFR) = (act(y(TGFbeta),rpar(:,125))*ymax(TGFR) - y(TGFR))/tau(TGFR); 
dydt(TIMP) = (OR(act(y(Akt),rpar(:,36)),OR(act(y(AT1R),rpar(:,50)),act(y(TGFR),rpar(:,128))))*ymax(TIMP) - y(TIMP))/tau(TIMP); 
dydt(TSC12) = (OR(inhib(y(Akt),rpar(:,12)),OR(inhib(y(ERK12),rpar(:,22)),OR(inhib(y(p38),rpar(:,26)),OR(inhib(y(RSK1),rpar(:,27)),OR(act(y(AMPK),rpar(:,38)),act(y(GSK3),rpar(:,69)))))))*ymax(TSC12) - y(TSC12))/tau(TSC12); 
dydt(Wnt) = (act(y(Wnt5),rpar(:,132))*ymax(Wnt) - y(Wnt))/tau(Wnt); 
dydt(Wnt5) = (OR(act(y(AT1R),rpar(:,51)),OR(act(y(NFkB),rpar(:,94)),OR(act(y(STAT),rpar(:,119)),act(y(TGFR),rpar(:,129)))))*ymax(Wnt5) - y(Wnt5))/tau(Wnt5); 

% utility functions 
function fact = act(x,rpar) 
% hill activation function with parameters w (weight), n (Hill coeff), EC50 
    w = rpar(1); 
    n = rpar(2); 
    EC50 = rpar(3); 
    beta = (EC50.^n - 1)./(2*EC50.^n - 1); 
    K = (beta - 1).^(1./n); 
    fact = w.*(beta.*x.^n)./(K.^n + x.^n); 
    if fact>w,                 % cap fact(x)<= 1 
        fact = w; 
    end
 
function finhib = inhib(x,rpar) 
% inverse hill function with parameters w (weight), n (Hill coeff), EC50 
    finhib = rpar(1) - act(x,rpar);
 
function z = OR(x,y) 
% OR logic gate 
    z = x + y - x*y;
 
function z = AND(rpar,varargin) 
% AND logic gate, multiplying all of the reactants together 
    w = rpar(1); 
    if w == 0, 
        z = 0; 
    else 
        v = cell2mat(varargin); 
        z = prod(v)/w^(nargin-2);  
    end 
