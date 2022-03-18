import math

r=1.85
L=41
rho_f=1800
rho_m=1170

def findMass(m):
	Mass=2*math.pi*r*L*m["H"]*(m["VR"]*rho_f+(1-m["VR"])*rho_m)*(2*m["h_r"]+m["c_r"]-m["h_r"]**2-2*m["h_r"]*m["c_r"])
	print(Mass)

modelA={
	"h_r":.120,
	"c_r":.132,
	"H":0.108,
	"VR":.875
}
	
modelB={
	"h_r" : .5,
	"c_r" : .265,
	"H" : 0.033,
	"VR" : .003
}

findMass(modelA)

