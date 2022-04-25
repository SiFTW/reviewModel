#file edited by reorderModel.jl to reorder the model species.

#######################################################
# Generated programmatically by CSV2JuliaDiffEq.      #
# http://github.com/SiFTW/CSV2JuliaDiffEq             #
#######################################################
# generated from:
#    reactions file: reactions.csv
#    parameters file file: parameters.csv
#    rate law file: rateLaws.csv
#
# Statistics:
#    Equations:26
#    Parameters:50
#######################################################



function odeFile(dy,y,p,t)
	LPS=maximum([y[1],0])
	CD14=maximum([y[2],0])
	CD14LPS=maximum([y[3],0])
	TLR4=maximum([y[4],0])
	TLR4LPS=maximum([y[5],0])
	CD14LPSen=maximum([y[6],0])
	TLR4en=maximum([y[7],0])
	TLR4LPSen=maximum([y[8],0])
	MYD88=maximum([y[9],0])
	MYD88s=maximum([y[10],0])
	TRIF=maximum([y[11],0])
	TRIFs=maximum([y[12],0])
	TRAF6=maximum([y[13],0])
	TRAF6s=maximum([y[14],0])
	IKKK_off=maximum([y[15],0])
	IKKK_on=maximum([y[16],0])
	IKK_off=maximum([y[17],0])
	IKK_on=maximum([y[18],0])
	IKK_i=maximum([y[19],0])
	TBK1=maximum([y[20],0])
	TBK1s=maximum([y[21],0])
	IRF3=maximum([y[22],0])
	IRF3s=maximum([y[23],0])
	IRF3n=maximum([y[24],0])
	IRF3ns=maximum([y[25],0])
	cbswitch=maximum([y[26],0])
	#LPS
	dy[1]= -CD14*LPS*8.75 + CD14LPS*0.07
	#CD14
	dy[2]= -CD14*LPS*8.75 + CD14LPS*0.07 + 0.00112 - CD14*0.00087793
	#CD14LPS
	dy[3]= + CD14*LPS*8.75 - CD14LPS*0.07 - TLR4*CD14LPS*5.543 + TLR4LPS*0.027715 - CD14LPS*0.065681 + CD14LPSen*0.04
	#TLR4
	dy[4]= -TLR4*CD14LPS*5.543 + TLR4LPS*0.027715 - TLR4*0.028028 + TLR4en*0.5 + 0.000052477
	#TLR4LPS
	dy[5]= + TLR4*CD14LPS*5.543 - TLR4LPS*0.027715 - TLR4LPS*0.065681 + TLR4LPSen*0.04
	#CD14LPSen
	dy[6]= -TLR4en*CD14LPSen*5.543 + TLR4LPSen*0.027715 + CD14LPS*0.065681 - CD14LPSen*0.04 - CD14LPSen*0.07
	#TLR4en
	dy[7]= -TLR4en*CD14LPSen*5.543 + TLR4LPSen*0.027715 + TLR4*0.028028 - TLR4en*0.5 - TLR4en*0.0653
	#TLR4LPSen
	dy[8]= + TLR4en*CD14LPSen*5.543 - TLR4LPSen*0.027715 + TLR4LPS*0.065681 - TLR4LPSen*0.04 - TLR4LPSen*0.0292 - TLR4LPSen*2*cbswitch
	#MYD88
	dy[9]= -(150*MYD88*(TLR4LPS^3))/((TLR4LPS^3)+(0.012448^3)) + MYD88s*2200 - MYD88*p
	#MYD88s
	dy[10]= + (150*MYD88*(TLR4LPS^3))/((TLR4LPS^3)+(0.012448^3)) - MYD88s*2200 + MYD88*p
	#TRIF
	dy[11]= -TRIF*18*TLR4LPSen + TRIFs*0.04
	#TRIFs
	dy[12]= + TRIF*18*TLR4LPSen - TRIFs*0.04
	#TRAF6
	dy[13]= -TRAF6*550*MYD88s - TRAF6*16*TRIFs + TRAF6s*18
	#TRAF6s
	dy[14]= + TRAF6*550*MYD88s + TRAF6*16*TRIFs - TRAF6s*18
	#IKKK_off
	dy[15]= -IKKK_off*0.098*TRAF6s + IKKK_on*2.5
	#IKKK_on
	dy[16]= + IKKK_off*0.098*TRAF6s - IKKK_on*2.5
	#IKK_off
	dy[17]= -IKK_off*1000*IKKK_on + IKK_on*0.9 + IKK_on*0 + IKK_i*2
	#IKK_on
	dy[18]= + IKK_off*1000*IKKK_on - IKK_on*0.9 - IKK_on*0
	#IKK_i
	dy[19]= -IKK_i*2
	#TBK1
	dy[20]= -TBK1*0.011*TRIFs + TBK1s*0.0358
	#TBK1s
	dy[21]= + TBK1*0.011*TRIFs - TBK1s*0.0358
	#IRF3
	dy[22]= -IRF3*25*TBK1s + IRF3s*0.010142 - IRF3*0.00091935 + IRF3n*3.6871 + 10.0546 - IRF3*0.0013
	#IRF3s
	dy[23]= + IRF3*25*TBK1s - IRF3s*0.010142 - IRF3s*0.0275 + IRF3ns*0.4145 - IRF3s*0.0584
	#IRF3n
	dy[24]= -IRF3n*0.0143 + IRF3ns*0.0194 + IRF3*0.00091935 - IRF3n*3.6871 - IRF3n*0.0091
	#IRF3ns
	dy[25]= + IRF3n*0.0143 - IRF3ns*0.0194 + IRF3s*0.0275 - IRF3ns*0.4145 - IRF3ns*0.015
	#cbswitch
	dy[26]=0
end
