#######################################################
# Generated programmatically by CSV2JuliaDiffEq.      #
# http://github.com/SiFTW/CSV2JuliaDiffEq             #
#######################################################
# generated from:
#    reactions file: reactions_myd88SELFACTIVATING.csv
#    parameters file file: parameters_myd88SELFACTIVATING.csv
#    rate law file: rateLaws.csv
#
# Statistics:
#    Equations:26
#    Parameters:49
#######################################################


using DifferentialEquations, Plots
pwd()

# function plotallspeciesn(phase2_sollist,phase1_sollist,tspan)
function plotallspeciesn(phase2_sollist,tspan)
	specie_names=["LPS","CD14","CD14LPS","TLR4pm","TLR4LPSpm","CD14LPSen","TLR4en","TLR4LPSen","MYD88","MYD88*","TRIF","TRIF*","TRAF6","TRAF6*","IKKK","IKKK*","IKK","IKK*","IKKi","TBK1","TBK1*","IRF3","IRF3*","IRF3n","IRF3n*"]
    subplotlist=[]
    if length(phase2_sollist) != length(phase1_sollist)
        error: ("unequal number of phase1 and phase2 solution vectors")
    end
    yvaluelistoflists=[]

    for j = 1:25
        yvaluelistoflists=[]
        for iter in range(1,length(phase2_sollist),step=1)
            yvalues=[]
            # for h in range(9900,10000,step=1) #100 timepoints pre stimulus are taken
            #     b=phase1_sollist[iter](h)
            #     push!(yvalues,b[j])
            # end
			if tspan==false
				for i in sol.t
					local a=phase2_sollist[iter](i)
					push!(yvalues,a[j])
				end
			else
	            for i in range(tspan[1],tspan[2],step=1)
	                local a=phase2_sollist[iter](i)
	                push!(yvalues,a[j])
	            end
			end
            push!(yvaluelistoflists,yvalues)
        end
         #make graph
		# x=tspan[1]:tspan[2]
		colorMap=Dict("WT"=>RGB(0,191/256,0),"L265P"=>RGB(128/256,0,128/256))
		colours=colorMap
		plot(yvaluelistoflists,title=specie_names[j],legend=:left,color_palette=[RGB(0,191/256,0),RGB(128/256,0,128/256)])
		# xlims!(tspan[1],tspan[2])
        push!(subplotlist,current()) #most recent graph appended
    end
    return subplotlist
end
function barallspeciesfinaln(phase2_sollist)
	specie_names=["LPS","CD14","CD14LPS","TLR4pm","TLR4LPSpm","CD14LPSen","TLR4en","TLR4LPSen","MYD88","MYD88*","TRIF","TRIF*","TRAF6","TRAF6*","IKKK","IKKK*","IKK","IKK*","IKKi","TBK1","TBK1*","IRF3","IRF3*","IRF3n","IRF3n*"]
    subplotlist=[]
    yvaluelistoflists=[]
	sims=["WT","L265P"]

	# colorMap=Dict("WT"=>RGB(0/256,100/256,0/256),"L265P"=>RGB(200/256,191/256,230/256))
	colorMap=Dict("WT"=>RGB(0,191/256,0),"L265P"=>RGB(128/256,0,128/256))
	colours=[colorMap[i] for i in sims]
    for j = 1:25
        finalvalues=[]
        for iter in range(1,length(phase2_sollist),step=1)
	        local a=phase2_sollist[iter][end][j]
            push!(finalvalues,a)
        end
		ub=maximum(finalvalues)*1.2
		lb=-0.2*ceil(ub,sigdigits=1)
		ylimz=(lb,ceil(ub,sigdigits=1))
		# ticks = round.(range(0,ub,length = 5),digits=1)
		ticks = (range(lb,ceil(ub,sigdigits=1),length = 3), string.(round.(range(lb,ceil(ub,sigdigits=1),length = 3),sigdigits=2)))
		ticks[2][end] = string(ceil(ub,sigdigits=1))
		if maximum(ticks[1])<1e-6
			lb=ceil(1e-10,sigdigits=1)
			ticks=[-2e-11, 0, 1e-10],["-2e-11","0.0","1e-10"]
			ylimz=(-2e-11,1e-10)
			finalvalues[end-1]=1e-13
			finalvalues[end]=1e-13
			bar(["WT","L265P"],finalvalues,colour=colours,title=specie_names[j],legend=:left,yticks=ticks,ylims=ylimz,ytickfontsize=15)
		else
		bar(["WT","L265P"],finalvalues,colour=colours,title=specie_names[j],legend=:left#=,yticks=ticks=#,ylims=ylimz,ytickfontsize=15)
		end
		# bar([finalvalues,title=specie_names[j],legend=:left)
        push!(subplotlist,current()) #most recent graph appended
    end
    return subplotlist
end

function justTLR_Myd88_SELFACTIVATING(dy,y,p,t)
	LPS=y[1]
	CD14=y[2]
	CD14LPS=y[3]
	TLR4pm=y[4]
	TLR4LPSpm=y[5]
	CD14LPSen=y[6]
	TLR4en=y[7]
	TLR4LPSen=y[8]
	MYD88=y[9]
	MYD88_A=y[10]
	TRIF=y[11]
	TRIF_A=y[12]
	TRAF6=y[13]
	TRAF6_A=y[14]
	IKKK=y[15]
	IKKK_A=y[16]
	IKK=y[17]
	IKK_A=y[18]
	IKKi=y[19]
	TBK1=y[20]
	TBK1_A=y[21]
	IRF3=y[22]
	IRF3_A=y[23]
	IRF3n=y[24]
	IRF3n_A=y[25]
	cbswitch=y[26]
	#LPS
	# dy[1]= -LPS*CD14*pmodfun(1) + CD14LPS*pmodfun(2)
	dy[1]= -0
	#CD14
	dy[2]= -LPS*CD14*pmodfun(1) + CD14LPS*pmodfun(2) + pmodfun(3) - CD14*pmodfun(4)
	#CD14LPS
	dy[3]= + LPS*CD14*pmodfun(1) - CD14LPS*pmodfun(2) - CD14LPS*TLR4pm*pmodfun(5) + TLR4LPSpm*pmodfun(6) - CD14LPS*pmodfun(10) + CD14LPSen*pmodfun(11)
	#TLR4pm
	dy[4]= -CD14LPS*TLR4pm*pmodfun(5) + TLR4LPSpm*pmodfun(6) + pmodfun(9) - TLR4pm*pmodfun(12) + TLR4en*pmodfun(13)
	#TLR4LPSpm
	dy[5]= + CD14LPS*TLR4pm*pmodfun(5) - TLR4LPSpm*pmodfun(6) - TLR4LPSpm*pmodfun(14) + TLR4LPSen*pmodfun(15)
	#CD14LPSen
	dy[6]= -CD14LPSen*TLR4en*pmodfun(7) + TLR4LPSen*pmodfun(8) + CD14LPS*pmodfun(10) - CD14LPSen*pmodfun(11) - CD14LPSen*pmodfun(16)
	#TLR4en
	dy[7]= -CD14LPSen*TLR4en*pmodfun(7) + TLR4LPSen*pmodfun(8) + TLR4pm*pmodfun(12) - TLR4en*pmodfun(13) - TLR4en*pmodfun(17)
	#TLR4LPSen
	dy[8]= + CD14LPSen*TLR4en*pmodfun(7) - TLR4LPSen*pmodfun(8) + TLR4LPSpm*pmodfun(14) - TLR4LPSen*pmodfun(15) - TLR4LPSen*pmodfun(18) - TLR4LPSen*pmodfun(19)*cbswitch
	#MYD88
	dy[9]= -(pmodfun(20)*MYD88*(TLR4LPSpm^pmodfun(21)))/((TLR4LPSpm^pmodfun(21))+pmodfun(22)^pmodfun(21)) + MYD88_A*pmodfun(23) - MYD88*pmodfun(50)
	#MYD88_A
	dy[10]= + (pmodfun(20)*MYD88*(TLR4LPSpm^pmodfun(21)))/((TLR4LPSpm^pmodfun(21))+(pmodfun(22)^pmodfun(21))) - MYD88_A*pmodfun(23) + MYD88*pmodfun(50)
	#TRIF
	dy[11]= -TRIF*TLR4LPSen*pmodfun(24) + TRIF_A*pmodfun(25)
	#TRIF_A
	dy[12]= + TRIF*TLR4LPSen*pmodfun(24) - TRIF_A*pmodfun(25)
	#TRAF6
	dy[13]= -TRAF6*MYD88_A*pmodfun(26) - TRAF6*TRIF_A*pmodfun(27) + TRAF6_A*pmodfun(28)
	#TRAF6_A
	dy[14]= + TRAF6*MYD88_A*pmodfun(26) + TRAF6*TRIF_A*pmodfun(27) - TRAF6_A*pmodfun(28)
	#IKKK
	dy[15]= -IKKK*TRAF6_A*pmodfun(29) + IKKK_A*pmodfun(30)
	#IKKK_A
	dy[16]= + IKKK*TRAF6_A*pmodfun(29) - IKKK_A*pmodfun(30)
	#I6]
	dy[17]= -IKK*IKKK_A*pmodfun(31) + IKK_A*pmodfun(32) + IKKi*pmodfun(34)
	#IKK_A
	dy[18]= + IKK*IKKK_A*pmodfun(31) - IKK_A*pmodfun(32) - IKK_A*pmodfun(33)
	#IKKi
	dy[19]= + IKK_A*pmodfun(33) - IKKi*pmodfun(34)
	#TBK1
	dy[20]= -TBK1*TRIF_A*pmodfun(35) + TBK1_A*pmodfun(36)
	#TBK1_A
	dy[21]= + TBK1*TRIF_A*pmodfun(35) - TBK1_A*pmodfun(36)
	#IRF3
	dy[22]= -IRF3*TBK1_A*pmodfun(37) + IRF3_A*pmodfun(38) - IRF3*pmodfun(41) + IRF3n*pmodfun(42) + pmodfun(45) - IRF3*pmodfun(46)
	#IRF3_A
	dy[23]= + IRF3*TBK1_A*pmodfun(37) - IRF3_A*pmodfun(38) + IRF3n_A*pmodfun(44) - IRF3_A*pmodfun(43) + IRF3n_A*pmodfun(44) - IRF3_A*pmodfun(47)
	#IRF3n
	dy[24]= -IRF3n*pmodfun(39) + IRF3*pmodfun(41) - IRF3n*pmodfun(42) - IRF3n*pmodfun(49)
	#IRF3n_A
	dy[25]= + IRF3n*pmodfun(39) - IRF3n_A*pmodfun(44) + IRF3_A*pmodfun(43) - IRF3n_A*pmodfun(44) - IRF3n_A*pmodfun(48)
	#cbswitch
	dy[26]= - 0
end
function condition(cheng2015_justTLR_variableparams,t,integrator)
   t-0
end
function affect!(integrator)
    #cbswitch for delayed TLR4LPSen degradation
    integrator.u[26]=1
end


u0 = [0, 0.08, 0, 0, 0.0001, 0, 0, 0, 0.1, 0, 0.1, 0, 0.1, 0, 0.1, 0, 0.1, 0, 0, 0.1, 0, 0.1, 0, 0, 0, 0]
u0[17]=0.1
u0test=zeros(26)
LPSinit=0.000
tspan=(0.0,60*18)

pfun=x->(8.75, 0.07, 0.00112, 0.000878, 5.543, 0.027715, 5.543, 0.027715, 5.25e-05, 0.065681, 0.04, 0.028028, 0.5, 0.065681, 0.04, 0.07, 0.0653, 0.0292, 2.0, 150.0, 3.0, 0.012448, 2200.0, 18.0, 0.04, 550.0, 16.0, 18.0, 0.098, 2.5, 1000.0, 0.9, 0.0, 2.0, 0.011, 0.0358, 25.0, 0.0101, 0.0143, 0.4145, 0.000919, 3.6871, 0.0275, 0.4145, 10.0546, 0.0013, 0.0584, 0.015, 0.0091, 1, 120)[x]
modifiers=ones(51)
pmodfun=x->(modifiers[x]*pfun(x))

sollist=[]
phase1_sollist=[]

for myd88_selfactivation in [0,500]
	modifiers[50]=myd88_selfactivation
	phase1_prob=DiffEqBase.ODEProblem(justTLR_Myd88_SELFACTIVATING,u0test,(0.0,100),pmodfun;maxiters=1e6)
	phase1_sol=DiffEqBase.solve(phase1_prob)
	# u0_phase2=phase1_sol[end]
	u0_phase2=u0
	u0_phase2[1]=LPSinit
	prob=DiffEqBase.ODEProblem(justTLR_Myd88_SELFACTIVATING, u0_phase2, tspan, pmodfun; maxiters=1e6)

	cb=ContinuousCallback(condition,affect!,nothing)
	cbs=CallbackSet(cb)
	sol=DiffEqBase.solve(prob, callback=cb)
	push!(phase1_sollist, phase1_sol)
	push!(sollist, sol)
end


subplotlist=plotallspeciesn(sollist,tspan)
plot(subplotlist[1],subplotlist[2],subplotlist[3],subplotlist[4],subplotlist[5],subplotlist[6],subplotlist[7],subplotlist[8],subplotlist[9],subplotlist[10],subplotlist[11],subplotlist[12],subplotlist[13],subplotlist[14],subplotlist[15],subplotlist[16],subplotlist[17],subplotlist[18],subplotlist[19],subplotlist[20],subplotlist[21],subplotlist[22],subplotlist[23],subplotlist[24],subplotlist[25],size=[1250,800],linewidth=1,show=true,legend=:none)
current()

testplots=barallspeciesfinaln(sollist)
plot(testplots[1],testplots[2],testplots[3],testplots[4],testplots[5],testplots[6],testplots[7],testplots[8],testplots[9],testplots[10],testplots[11],testplots[12],testplots[13],testplots[14],testplots[15],testplots[16],testplots[17],testplots[18],testplots[19],testplots[20],testplots[21],testplots[22],testplots[23],testplots[24],testplots[25],size=[1250,1800]*0.6,legend=:none)
