#! /usr/bin/python
import os

#IOV_list= ['2016P8','2016QCD','2016BCD','2016EF',
#           '2016APVP8','2016APVQCD','2016FGH',
#           '2017P8','2017QCD','2017B','2017C','2017D','2017E','2017F',
#           '2018P8','2018QCD','2018A1','2018A2','2018B','2018C',
#           '2018D1','2018D2','2018D3','2018D4']
IOV_list= [
    'UL2016BCD','UL2016EF','UL2016GH',
    'UL2017B','UL2017C','UL2017D','UL2017E','UL2017F',
    'UL2018A','UL2018B','UL2018C',#'UL2018D',
    'UL2018D1','UL2018D2',
    'UL2016APVMG','UL2016MG','UL2018MG', #'UL2017MG'
]
version = 'v25'

#os.system("rm *.so *.d *.pcm")
os.system("root -l -b -q mk_CondFormats.C")
for iov in IOV_list:
    print "Process DijetHistosFill.C+g for IOV "+iov
    os.system("ls -ltrh rootfiles/jmenano_mc_out_"+iov+"_"+version+".root")
    os.system("ls -ltrh rootfiles/jmenano_data_out_"+iov+"_"+version+".root")
    os.system("ls -ltrh logs/log_"+iov+"_"+version+".txt")
    os.system("nohup root -l -b -q 'mk_DijetHistosFill.C(\""+iov+"\",\""+version+"\")' > logs/log_"+iov+"_"+version+".txt &")
    print " => Follow logging with 'tail -f logs/log_"+iov+"_"+version+".txt'"
#    os.system("fs flush")
#    wait()
#    time.sleep(sleep_time)
