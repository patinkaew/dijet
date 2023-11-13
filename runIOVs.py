#! /usr/bin/python
import os
import argparse

IOV_list= [    'UL2016BCD','UL2016EF','UL2016GH',
   'UL2017B','UL2017C','UL2017D','UL2017E','UL2017F',
    'UL2018A','UL2018B','UL2018C', 'UL2018D',
    'UL2018D1','UL2018D2',
    'UL2016APVMG','UL2016MG','UL2017MG', 'UL2018MG',
    'UL2016BCD_ZB','UL2016EF_ZB','UL2016GH_ZB',
    'UL2017B_ZB','UL2017C_ZB','UL2017D_ZB','UL2017E_ZB','UL2017F_ZB',
    'UL2018A_ZB','UL2018B_ZB','UL2018C_ZB', 'UL2018D_ZB',
    '2022C','2022D','2022E','2022F1','2022F2','2022G',
    '2022C_ZB','2022D_ZB','2022E_ZB','2022F_ZB','2022G_ZB',
    '2023BCv123','2023Cv4','2023D',
    '2023BCv123_ZB','2023Cv4_ZB','2023D_ZB',
    'Summer22MG','Summer22EEMG',
    'Summer22MG1','Summer22MG2',
    'Summer22EEMG1','Summer22EEMG2','Summer22EEMG3','Summer22EEMG4'
]
version = 'v32'

IOV_input = []

parser = argparse.ArgumentParser(description='Run all IOVs')

# The user can pass the IOV list, version and max number of files as an argument
parser.add_argument('--IOV_list', nargs='+', default=IOV_input)
parser.add_argument('--version', default=version)
parser.add_argument('--max_files', default=9999)
args = parser.parse_args()

if args.IOV_list:
    # if the user passes all, then all IOVs are run
    if 'all' in args.IOV_list:
        IOV_input = IOV_list
    else:
        # Check that all IOVs passed are in the list
        for iov in args.IOV_list:
            if iov not in IOV_list:
                print('IOV '+iov+' not in list of IOVs')
                exit()
            else:
                IOV_input.append(iov)
else:
    print('No IOV list passed')
    exit()
    
print('IOVs to run: ', IOV_input)

#os.system("rm *.so *.d *.pcm")
os.system("root -l -b -q mk_CondFormats.C")
for iov in IOV_input:
    print(f"Process DijetHistosFill.C+g for IOV {iov}")
    os.system(f"ls -ltrh rootfiles/jmenano_mc_out_{iov}_{version}.root")
    os.system(f"ls -ltrh rootfiles/jmenano_data_out_{iov}_{version}.root")
    os.system(f"ls -ltrh logs/log_{iov}_{version}.txt")
    os.system(f"nohup root -l -b -q 'make/mk_DijetHistosFill.C(\"{iov}\",\"{version}\",{args.max_files})' > logs/log_{iov}_{version}.txt &")
    print(f" => Follow logging with 'tail -f logs/log_{iov}_{version}.txt'")
#    os.system("fs flush")
#    wait()
#    time.sleep(sleep_time)
