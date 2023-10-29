#! /usr/bin/python
import os

#IOV_list= ['2016P8','2016QCD','2016BCD','2016EF',
#           '2016APVP8','2016APVQCD','2016FGH',
#           '2017P8','2017QCD','2017B','2017C','2017D','2017E','2017F',
#           '2018P8','2018QCD','2018A1','2018A2','2018B','2018C',
#           '2018D1','2018D2','2018D3','2018D4']
IOV_list= [
#    'UL2016BCD','UL2016EF','UL2016GH',
#    'UL2017B','UL2017C','UL2017D','UL2017E','UL2017F',
#    'UL2018A','UL2018B','UL2018C',#'UL2018D',
#    'UL2018D1','UL2018D2',
#    'UL2016APVMG','UL2016MG','UL2017MG', 'UL2018MG'
#    'UL2016BCD_ZB','UL2016EF_ZB','UL2016GH_ZB',
#    'UL2017B_ZB','UL2017C_ZB','UL2017D_ZB','UL2017E_ZB','UL2017F_ZB',
#    'UL2018A_ZB','UL2018B_ZB','UL2018C_ZB', 'UL2018D_ZB'
    '2022C','2022D','2022E','2022F1','2022F2','2022G',
    '2022C_ZB','2022D_ZB','2022E_ZB','2022F_ZB','2022G_ZB',
    '2023BCv123','2023Cv4','2023D',
    '2023BCv123_ZB','2023Cv4_ZB','2023D_ZB',
#    'Summer22MG','Summer22EEMG'
    'Summer22MG1','Summer22MG2',
    'Summer22EEMG1','Summer22EEMG2','Summer22EEMG3','Summer22EEMG4'
]
version = 'v32'

#os.system("rm *.so *.d *.pcm")
os.system("root -l -b -q mk_CondFormats.C")
for iov in IOV_list:
    print("Process DijetHistosFill.C+g for IOV "+iov)
    os.system("ls -ltrh rootfiles/jmenano_mc_out_"+iov+"_"+version+".root")
    os.system("ls -ltrh rootfiles/jmenano_data_out_"+iov+"_"+version+".root")
    os.system("ls -ltrh logs/log_"+iov+"_"+version+".txt")
    os.system("nohup root -l -b -q 'mk_DijetHistosFill.C(\""+iov+"\",\""+version+"\")' > logs/log_"+iov+"_"+version+".txt &")
    print(" => Follow logging with 'tail -f logs/log_"+iov+"_"+version+".txt'")
#    os.system("fs flush")
#    wait()
#    time.sleep(sleep_time)
