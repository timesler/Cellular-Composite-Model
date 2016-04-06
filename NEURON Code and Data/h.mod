TITLE I-h channel from Magee 1998

NEURON {
	SUFFIX hd
	NONSPECIFIC_CURRENT i
    RANGE ghdbar, vhalfl, ghd, i, kl
    GLOBAL qtl
}
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 				(mV)
    ehd  			(mV)
	celsius 		(degC)
	ghdbar=1e-7	(mho/cm2)
    vhalfl=-81   	(mV)
	kl=-8
	vhalft=-75   	(mV)
	a0t=0.011      	(/ms)
	zetat=2.2    	(1)
	gmt=.4   		(1)
	q10=4.5
	qtl=1
}

STATE {
        l
}

ASSIGNED {
	i (mA/cm2)
	linf
	taul
	ghd
}

INITIAL {
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ghd = ghdbar*l
	i = ghd*(v-ehd)
}

DERIVATIVE states {
        rate(v)
        l' =  (linf - l)/taul
}

FUNCTION alpt(v(mV)) {
  		alpt = exp(0.0378*zetat*(v-vhalft))
}

FUNCTION bett(v(mV)) {
  		bett = exp(0.0378*zetat*gmt*(v-vhalft))
}


PROCEDURE rate(v (mV)) {
        LOCAL a,qt
        qt=q10^((celsius-33)/10)
        a = alpt(v)
        linf = 1/(1 + exp(-(v-vhalfl)/kl))
        taul = bett(v)/(qtl*qt*a0t*(1+a))
}














