import numpy as np
import pandas as pd
import os
import math 
from itertools import combinations


def select_genotypes(df, sample_list):
    filter_df = df.set_index(['seqname', 'start', 'end']).filter(items=sample_list, axis=1)
    return filter_df.reset_index()
def file_to_list(file):
    file_l = pd.read_csv(file)['genotype'].tolist()
    return file_l



#for each genome the script take aroun 12 minutes
for namels in 'landmark','chinese_spring','CWI86942':
    
    #importing the data to support the Aegilops tauschii population structure
    in_file = pd.read_csv(f'./data/{namels}_995_combined_queries_50000w.tsv.gz', delimiter='\t')
    my_list_231 = file_to_list('./genotypes_K22.tsv')
    my_list_3 = file_to_list('./genotypes_L3.tsv')
    all_lin231 = select_genotypes(in_file, my_list_231)
    all_lin3 = select_genotypes(in_file, my_list_3)
    pops_L231 = pd.read_csv('./allchrv14_pops.csv', delimiter=',', index_col=0, header=None)
    alll_llin231 = all_lin231[all_lin231.columns[3:]]
    alll_llin3 = all_lin3[all_lin3.columns[3:]]
    
    
    my_list_1 = file_to_list('./genotypes_L1_K22.tsv')
    my_list_2 = file_to_list('./genotypes_L2_K22.tsv')
    all_lin1 = select_genotypes(in_file, my_list_1)
    all_lin2 = select_genotypes(in_file, my_list_2)
    alll_llin1 = all_lin1[all_lin1.columns[3:]]
    alll_llin2 = all_lin2[all_lin2.columns[3:]]


    predf_ff= {'chrom':[],'start':[],'lineage':[],'K9_population':[],'lmi':[],'lmit':[],'sll':[],'soccf':[],'K15_population':[],'soccf_15':[],'K22_population':[],'soccf_22':[]}
    resu_Df_ff = pd.DataFrame(predf_ff)

    ltot =len(all_lin231)

    for i in range(ltot):
        datad = []
        diction = {}
        ollistlat = []
        ollistlut = []
        ollistlot = []
        spop = ''
        soccf = 0
        sspop = ''
        sll = 0
        ssspop = ''
        sconf = 1

		#take the position
        first_half = all_lin231.iloc[i, :2].values.flatten().tolist()
        #take the names of the accessions with IBSpy variations score lower than 31
        idbs_acc231 = alll_llin231.columns[alll_llin231.loc[i].lt(31)].tolist()
        #names of accessions to determine undetermined lineage 3
        idbs_acc3 = alll_llin3.columns[alll_llin3.loc[i].lt(216)].tolist()



        if idbs_acc231:
            linnn231 = True
        else:
            linnn231 = False
        if idbs_acc3:
            linnn3 = True
        else:
            linnn3 = False

        if linnn231 == False and linnn3 == False:
        	#determine if the window is an undetermined lineage 1 or 2 or a non tauschii derived region
            idbs_acc2_br_g = alll_llin2.columns[alll_llin2.loc[i].lt(256)].tolist()
            idbs_acc1_br_g = alll_llin1.columns[alll_llin1.loc[i].lt(256)].tolist()
            idbs_acc2_br_p = alll_llin2.columns[alll_llin2.loc[i].gt(30)].tolist()
            idbs_acc1_br_p = alll_llin1.columns[alll_llin1.loc[i].gt(30)].tolist()
            idbs_acc2_br = len(list(set(idbs_acc2_br_g).intersection(idbs_acc2_br_p)))/271
            idbs_acc1_br = len(list(set(idbs_acc1_br_g).intersection(idbs_acc1_br_p)))/576

            if ((idbs_acc2_br >= idbs_acc1_br) and (idbs_acc2_br != 0)):
                linnn = 'L2'
                spop = 'UND_L2'
                soccf = 0
                soccf_15 = 0
                soccf_22 = 0
                lmi = 0
                lmit = 0
                sspop = 'UND_L2'
                ssspop = 'UND_L2'
            elif ((idbs_acc1_br > idbs_acc2_br) and (idbs_acc1_br != 0)):
                linnn = 'L1'
                spop = 'UND_L1'
                soccf = 0
                soccf_15 = 0
                soccf_22 = 0
                lmi = 0
                lmit = 0
                sspop = 'UND_L1'
                ssspop = 'UND_L1'
            else:
                linnn = 'NO_TAU'
                spop = 'NO_TAU'
                soccf = 0
                soccf_15 = 0
                soccf_22 = 0
                lmi = 0
                lmit = 0
                sspop = 'NO_TAU'
                ssspop = 'NO_TAU'  
        elif linnn231 == False and linnn3 == True:
            linnn = 'L3'
            spop = 'L3'
            #sconf = 1
            soccf = 0
            soccf_15 = 0
            soccf_22 = 0
            lmi = 0
            lmit = 0
            sspop = 'L3'
            ssspop = 'L3'
        elif linnn231 == True:
            linnn = 'L1/2'
        else:
            linnn = 'UNK'
            spop = 'UNK'
            #sconf = 1
            soccf = 0
            soccf_15 = 0
            soccf_22 = 0
            lmi = 0
            lmit = 0
            sspop = 'UNK'
            ssspop = 'UNK'




        #subpopulation determination for K=9
        flag = 0
        sll = 31
        #counting the number of identical-by-state accessions
        lmi = len(alll_llin231.columns[alll_llin231.loc[i].lt(31)].tolist())
        lmit = lmi
        for ll in 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31:
            if flag == 0:
                idbs_acc231 = alll_llin231.columns[alll_llin231.loc[i].lt(ll)].tolist()
                #using a root function as condition to pick an appropriate amount of accessions as a limit
                if len(idbs_acc231) > ((lmi**(1/2))*2)+10:
                    flag = 1
                    #upper boundary of IBSpy score used
                    sll = ll
                    #counting how many accessions are selected
                    lmit = len(idbs_acc231)

        for x in idbs_acc231:
            ollat = pops_L231.loc[x, 1]
            ollistlat.append(ollat)
            ollistlat[:] = (value for value in ollistlat if value != 'ADMIXED')
        if ollistlat:
            flag = 1

            arms = ollistlat.count('L1W')
            chis = ollistlat.count('L1E-Y')
            tjks = ollistlat.count('L1E-X')
            gils = ollistlat.count('L2E-1')
            caus = ollistlat.count('L2W-1')
            tkms = ollistlat.count('L2W-2')
            gors = ollistlat.count('L2E-2')
            li3s = ollistlat.count('L3')
            
            #tangent function computation to take in account variability of accessions number per subpopulation
            sgil = 2 * math.tan(((1/(54*2.3/math.pi*2))*(gils))+math.pi/1.7)
            scau = 2 * math.tan(((1/(111*2.3/math.pi*2))*(caus))+math.pi/1.7)
            stkm = 2 * math.tan(((1/(31*2.3/math.pi*2))*(tkms))+math.pi/1.7)
            sgor = 2 * math.tan(((1/(75*2.3/math.pi*2))*(gors))+math.pi/1.7)
            sarm = 2 * math.tan(((1/(97*2.3/math.pi*2))*(arms))+math.pi/1.7)
            schi = 2 * math.tan(((1/(71*2.3/math.pi*2))*(chis))+math.pi/1.7)
            stjk = 2 * math.tan(((1/(408*2.3/math.pi*2))*(tjks))+math.pi/1.7)
            sli3 = 2 * math.tan(((1/(11*2.3/math.pi*2))*(li3s))+math.pi/1.7)

            diction = {'L3': sli3, 'L2E-1': sgil, 'L2W-1': scau, 'L2W-2': stkm, 'L2E-1' : sgor, 'L1W': sarm, 'L1E-Y': schi, 'L1E-X' : stjk}
            slist_9 = [arms, chis, tjks, gils, caus, tkms, gors]
            #computing the sum of the pairwise products of the element of the list as a score for accuracy of the prediction, if only one population is present the score is 0
            soccf = sum(a*b for a,b in combinations(slist_9,2))
            #set a maximum score value to have nice plottings
            if soccf > 500:
                soccf = 500
            spop = max(diction, key=diction.get)
            #part of the script supporting the results in the paper ends here
            
            #the following part is predicting the subpopulation contributions based on higher number of K from ancestry analysis (K=15, K=22)
            #names of the subpopulations matches the geographical names in which most of the accessions from the subpopulation were collected
            
            for y in idbs_acc231:
                ollet = pops_L231.loc[y, 1]
                if ((ollet == spop) or (ollet == 'ADMIXED')):
                    ollot = pops_L231.loc[y, 2]
                    ollistlot.append(ollot)
            ollistlot[:] = (value for value in ollistlot if value != 'ADMIXED')

            if ollistlot:


                gilss = ollistlot.count('GILAN')
                causs = ollistlot.count('CAUCASUS')
                golss = ollistlot.count('GOLESTAN')
                gorss = ollistlot.count('GORGAN')
                mazss = ollistlot.count('MAZANDARAN')
                tkmss = ollistlot.count('TURKMENISTAN')
                urmss = ollistlot.count('URMIA')
                armss = ollistlot.count('ARMENIA')
                chnss = ollistlot.count('CHINA')
                geoss = ollistlot.count('GEORGIA')
                pakss = ollistlot.count('PAKISTAN')
                tajss = ollistlot.count('TAJIKISTAN')
                xinss = ollistlot.count('XINJIANG')
                li3ss = ollistlot.count('L3')
                
                ssarm = 2 * math.tan(((1/(42*2.3/math.pi*2))*(armss))+math.pi/1.7)
                sschn = 2 * math.tan(((1/(69*2.3/math.pi*2))*(chnss))+math.pi/1.7)
                ssgeo = 2 * math.tan(((1/(31*2.3/math.pi*2))*(geoss))+math.pi/1.7)
                sspak = 2 * math.tan(((1/(19*2.3/math.pi*2))*(pakss))+math.pi/1.7)
                sstaj = 2 * math.tan(((1/(76*2.3/math.pi*2))*(tajss))+math.pi/1.7)
                ssxin = 2 * math.tan(((1/(25*2.3/math.pi*2))*(xinss))+math.pi/1.7)
                ssgil = 2 * math.tan(((1/(37*2.3/math.pi*2))*(gilss))+math.pi/1.7)
                sscau = 2 * math.tan(((1/(110*2.3/math.pi*2))*(causs))+math.pi/1.7)
                ssgol = 2 * math.tan(((1/(36*2.3/math.pi*2))*(golss))+math.pi/1.7)
                ssgor = 2 * math.tan(((1/(6*2.3/math.pi*2))*(gorss))+math.pi/1.7)
                ssmaz = 2 * math.tan(((1/(16*2.3/math.pi*2))*(mazss))+math.pi/1.7)
                sstkm = 2 * math.tan(((1/(6*2.3/math.pi*2))*(tkmss))+math.pi/1.7)
                ssurm = 2 * math.tan(((1/(18*2.3/math.pi*2))*(urmss))+math.pi/1.7)
                ssli3 = 2 * math.tan(((1/(11*2.3/math.pi*2))*(li3ss))+math.pi/1.7)

                diction = {'L3': ssli3, 'ARMENIA': ssarm, 'CHINA': sschn, 'GEORGIA' : ssgeo, 'PAKISTAN' : sspak, 'TAJIKISTAN' : sstaj, 'XINJIANG' : ssxin, 'GILAN': ssgil, 'CAUCASUS': sscau, 'GOLESTAN' : ssgol, 'GORGAN' : ssgor, 'MAZANDARAN': ssmaz, 'TURKMENISTAN' : sstkm, 'URMIA' : ssurm}
                sspop = max(diction, key=diction.get)
                slist_15 = [gilss, causs, golss, gorss, mazss, tkmss, urmss, armss, chnss, geoss, pakss, tajss, xinss]
                soccf_15 = sum(a*b for a,b in combinations(slist_15,2))

                for z in idbs_acc231:
                    ollut = pops_L231.loc[z, 1]
                    if ((ollut == spop) or (ollut == 'ADMIXED')):
                        ollit = pops_L231.loc[z, 3]
                        ollistlut.append(ollit)
                ollistlut[:] = (value for value in ollistlut if value != 'ADMIXED')
                #print(ollistlut)
                if ollistlut:
                    

                    causss = ollistlut.count('CAUCASUS')
                    gilsss = ollistlut.count('GILAN')
                    goasss = ollistlut.count('GOLESTAN_A')
                    gobsss = ollistlut.count('GOLESTAN_B')
                    gorsss = ollistlut.count('GORGAN')
                    lansss = ollistlut.count('LANKARAN')
                    mazsss = ollistlut.count('MAZANDARAN')
                    mugsss = ollistlut.count('MUGHAN')
                    tkmsss = ollistlut.count('TURKMENISTAN')
                    urmsss = ollistlut.count('URMIA')
                    afgsss = ollistlut.count('AFGHANISTAN')
                    armsss = ollistlut.count('ARMENIA')
                    ascsss = ollistlut.count('ASIAN_MOUNTAIN_CORRIDOR')
                    ceasss = ollistlut.count('CENTRAL_ASIA')
                    chnsss = ollistlut.count('CHINA')
                    geosss = ollistlut.count('GEORGIA')
                    paksss = ollistlut.count('PAKISTAN')
                    tajsss = ollistlut.count('TAJIKISTAN')
                    trnsss = ollistlut.count('TURAN')
                    xinsss = ollistlut.count('XINJIANG')
                    li3sss = ollistlut.count('L3')

                    ssscau = 2 * math.tan(((1/(61*2.3/math.pi*2))*(causss))+math.pi/1.7)
                    sssgil = 2 * math.tan(((1/(12*2.3/math.pi*2))*(gilsss))+math.pi/1.7)
                    sssgoa = 2 * math.tan(((1/(16*2.3/math.pi*2))*(goasss))+math.pi/1.7)
                    sssgob = 2 * math.tan(((1/(16*2.3/math.pi*2))*(gobsss))+math.pi/1.7)
                    sssgor = 2 * math.tan(((1/(6*2.3/math.pi*2))*(gorsss))+math.pi/1.7)
                    ssslan = 2 * math.tan(((1/(16*2.3/math.pi*2))*(lansss))+math.pi/1.7)
                    sssmaz = 2 * math.tan(((1/(14*2.3/math.pi*2))*(mazsss))+math.pi/1.7)
                    sssmug = 2 * math.tan(((1/(9*2.3/math.pi*2))*(mugsss))+math.pi/1.7)
                    ssstkm = 2 * math.tan(((1/(6*2.3/math.pi*2))*(tkmsss))+math.pi/1.7)
                    sssurm = 2 * math.tan(((1/(18*2.3/math.pi*2))*(urmsss))+math.pi/1.7)
                    sssafg = 2 * math.tan(((1/(19*2.3/math.pi*2))*(afgsss))+math.pi/1.7)
                    sssarm = 2 * math.tan(((1/(41*2.3/math.pi*2))*(armsss))+math.pi/1.7)
                    sssasc = 2 * math.tan(((1/(8*2.3/math.pi*2))*(ascsss))+math.pi/1.7)
                    ssscea = 2 * math.tan(((1/(10*2.3/math.pi*2))*(ceasss))+math.pi/1.7)
                    ssschn = 2 * math.tan(((1/(69*2.3/math.pi*2))*(chnsss))+math.pi/1.7)
                    sssgeo = 2 * math.tan(((1/(30*2.3/math.pi*2))*(geosss))+math.pi/1.7)
                    ssspak = 2 * math.tan(((1/(12*2.3/math.pi*2))*(paksss))+math.pi/1.7)
                    ssstaj = 2 * math.tan(((1/(35*2.3/math.pi*2))*(tajsss))+math.pi/1.7)
                    ssstrn = 2 * math.tan(((1/(6*2.3/math.pi*2))*(trnsss))+math.pi/1.7)
                    sssxin = 2 * math.tan(((1/(22*2.3/math.pi*2))*(xinsss))+math.pi/1.7)
                    sssli3 = 2 * math.tan(((1/(22*2.3/math.pi*2))*(li3sss))+math.pi/1.7)

                    slist_22 = [causss, gilsss, goasss, gobsss, gorsss, lansss, mazsss, mugsss, tkmsss, urmsss, afgsss, armsss, ascsss, ceasss, chnsss, geosss, paksss, tajsss, trnsss, xinsss]
                    soccf_22 = sum(a*b for a,b in combinations(slist_22,2))
                    diction = {'L3': sssli3, 'AFGHANISTAN': sssafg, 'ARMENIA' : sssarm, 'ASIAN_MOUNTAIN_CORRIDOR' : sssasc, 'CENTRAL_ASIA' : ssscea, 'CHINA' : ssschn, 'GEORGIA' : sssgeo, 'PAKISTAN' : ssspak, 'TAJIKISTAN' : ssstaj, 'TURAN' : ssstrn, 'XINJIANG' : sssxin, 'CAUCASUS': ssscau, 'GILAN': sssgil, 'GOLESTAN_A' : sssgoa, 'GOLESTAN_B' : sssgob, 'GORGAN' : sssgor, 'LANKARAN': ssslan, 'MAZANDARAN': sssmaz, 'MUGHAN' : sssmug, 'TURKMENISTAN' : ssstkm, 'URMIA' : sssurm}
                    ssspop = max(diction, key=diction.get)


                else:
                    ssspop = 'UND_' + sspop
                    soccf_22 = 0
            else:
                sspop = 'UND_' + spop
                ssspop = 'UND_' + spop
                soccf_15 = 0
                soccf_22 = 0

        else:
            if linnn == 'L1':
                spop = 'UND_L1'
                sspop = 'UND_L1'
                ssspop = 'UND_L1'
                soccf = 0
                soccf_15 = 0
                soccf_22 = 0
            elif linnn == 'L2':
                spop = 'UND_L2'
                sspop = 'UND_L2'
                ssspop = 'UND_L2'
                soccf = 0
                soccf_15 = 0
                soccf_22 = 0
            elif linnn == 'L3':
                spop = 'UND_L3'
                sspop = 'UND_L3'
                ssspop = 'UND_L3'
                soccf = 0
                soccf_15 = 0
                soccf_22 = 0
            elif linnn == 'NO_TAU':
                spop = 'NO_TAU'
                sspop = 'NO_TAU'
                ssspop = 'NO_TAU'
                soccf = 0
                soccf_15 = 0
                soccf_22 = 0
            else:
                spop = 'UND'
                sspop = 'UND'
                ssspop = 'UND'
                soccf = 0
                soccf_15 = 0
                soccf_22 = 0




        datad = [linnn, spop, lmi, lmit, sll, soccf, sspop, soccf_15, ssspop, soccf_22]
        ddatad = first_half + datad

        resu_Df_ff.loc[len(resu_Df_ff)] = ddatad

    resu_Df_ff.to_csv(f'./{namels}_n231_5dk22.csv')

