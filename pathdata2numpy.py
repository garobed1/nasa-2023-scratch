"""
Convert file I/O path data from .att files to numpy arrays

Optionally save as npx, but can use internally as well

"""
import os, sys
import numpy as np
import json

def read_paths(source_dir, ext='.att', prefix='EDW_ERA5_', 
                props=['TEMP', 'WINDX', 'WINDY', 'HUMIDITY', 'PRESSURE'],  
                times=['01','02','03','04','05','06','07','08','09','10','11','12','13',
            '14','15','16','17','18','19','20','21','22','23'],
                days=['01','02','03','04','05','06','07','08','09','10','11','12','13',
            '14','15','16','17','18','19','20','21','22','23','24','25','26','27',
            '28','29','30','31'],
                months=['01','02','03','04','05','06','07','08','09','10','11', '12'],
                years=np.arange(2001,2023),
                filename='fulldata.json',
                write_to_json=False):

    # get names of files corresponding to dates and times
    valid_files = []
    for y in years:
        for m in months:
            for d in days:
                for h in times:
                    s = f'{prefix}{y}{m}{d}{h}{ext}'
                    valid_files.append(s)

    # count number of files
    attcounter = 0
    # list file names in chronological order, note first day
    daystrlist = []
    daynumlist = []
    for file in os.listdir(source_dir):
        if file in valid_files:
            attcounter += 1
            daystrlist.append(file)
            date = int(file.split('_')[-1].split('.')[0])
            daynumlist.append(date)
    print(f'{attcounter} paths as data')
    
    # reorder by date, then index numpy array based on list order
    newind = np.argsort(daynumlist)
    # print(newind)
    daynumlist[:] = [daynumlist[i] for i in newind]
    daystrlist[:] = [daystrlist[i] for i in newind]
    # print(daynumlist)
    # print(daystrlist)
    # return 0

    fulldata = {}
    fulldata['n'] = attcounter
    fulldata['filenames'] = daystrlist
    fulldata['dates'] = daynumlist
    fulldata['altitude'] = np.zeros([attcounter, 1])
    for name in props:
        fulldata[name] = np.zeros([attcounter, 1])

    # begin looping through files
    gfc = 0
    for attfile in daystrlist:
        
        with open(source_dir + '/' + attfile, 'r') as atf:
            # useful flags
            glc = 0
            lc = 0
            newvar = False
            propname = None

            # begin line loop
            for line in atf:
                if glc < 2: # skip first two lines
                    pass
                elif newvar: # initialize arrays, the line should be a single int
                    ndat = int(line)
                    diff = ndat - fulldata[propname].shape[1]
                    if diff > 0:
                        fulldata[propname] = np.pad(fulldata[propname], [(0, 0), (0, diff)], mode='constant')
                    diffa = ndat - fulldata['altitude'].shape[1]
                    if diffa > 0:
                        fulldata['altitude'] = np.pad(fulldata['altitude'], [(0, 0), (0, diff)], mode='constant')
                    newvar = False
                    lc = 0
                elif line[:1].isalpha(): # new var spotted
                    propname = line.rstrip()
                    if propname not in props:
                        Exception(f'Key {propname} does not match anything in props arg: {props}')
                    newvar = True
                else: # add data
                    dat = line.rstrip().split('\t')
                    fulldata['altitude'][gfc, lc] = float(dat[0])
                    fulldata[propname][gfc, lc] = float(dat[1])
                    lc += 1
                glc += 1
        gfc += 1

    if write_to_json:
        fulldata_list = {}
        for k,v in fulldata.items():
            if isinstance(v, np.ndarray):
                fulldata_list[k] = v.tolist() 
            else:
                fulldata_list[k] = v
        import pdb; pdb.set_trace()
        with open(filename, 'w') as fj:
            json.dump(fulldata_list, fj)

        return 0

    # find out if this is too memory intensive
    return fulldata


# test out function here
if __name__ == '__main__':
    source_dir = sys.argv[1]
    months = ['08']
    times = ['18']

    read_paths(source_dir, 
        months=months,
        times=times,
        filename=f'partial_{months[0]}{times[0]}_data.json',
        write_to_json=True)

     