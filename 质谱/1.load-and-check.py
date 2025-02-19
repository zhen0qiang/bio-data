from pyteomics import mgf, mass
import numpy as np
import math


cr = {1: 1, 2: 0.9, 3: 0.85, 4: 0.8, 5: 0.75, 6: 0.75, 7: 0.75, 8: 0.75}
types = {'un': 0, 'cid': 1, 'etd': 2, 'hcd': 3, 'ethcd': 4, 'etcid': 5}

def readmgf(fn, type):
    file = open(fn, "r")
    data = mgf.read(file, convert_arrays=1, read_charges=False, dtype='float32', use_index=False)

    # codes : list[dict] {'pep': pep-str, 'charge': c-int, 'mass': mass-float,
    # 'mz': mz-array[float], 'it': it-array[float], 'nmod': nmod-int, 'mod': mod-array[int]}
    codes = convert_mgf(data)

    for sp in codes:
        sp['type'] = types[type]
    return codes

def filter_spectra(db):
    return [sp for sp in db if spectra_ok(sp)]

def i2l(sps):
    sps = [sp.copy() for sp in sps]
    for sp in sps:
        sp['pep'] = sp['pep'].replace('I', 'L')
    return sps


def convert_mgf(sps):
    db = []

    for sp in sps:
        param = sp['params']

        if not 'charge' in param: raise
            
        c = int(str(param['charge'][0])[0])

        pep = title = param['title']
        if 'seq' in param: pep = param['seq']

        if 'pepmass' in param: mass = param['pepmass'][0]
        else: mass = float(param['parent'])

        rtime = 0 if not 'RTINSECONDS' in param else float(param['RTINSECONDS'])

        if 'hcd' in param:
            try:
                hcd = param['hcd']
                if hcd[-1] == '%':
                    hcd = float(hcd)
                elif hcd[-2:] == 'eV':
                    hcd = float(hcd[:-2])
                    hcd = hcd * 500 * cr[c] / mass
                else:
                    raise Exception("Invalid type!")
            except:
                hcd = 0
        else: hcd = 0

        mz = sp['m/z array']
        it = sp['intensity array']

        db.append({'pep': pep, 'charge':c, 'mass': mass, 'mz': mz, 'it': it, 'nmod': 0,
                   'mod': np.zeros(len(pep), 'int32'), # currently no mod supported
                   'nce': hcd, 'title': title })

    return db

def spectra_ok(sp, ppm_threshold=10): # check 
    mz, mass, pep, c = sp['mz'], sp['mass'], sp['pep'], sp['charge']

    if not pep.isalpha():
        return False # unknown mod

    if ppm_threshold > 0 and abs(ppmdiff(sp)) > ppm_threshold:
        return False

    return True

# 定义一个函数ppmdiff，用于计算质谱中肽的质量与理论质量的偏差
def ppmdiff(sp, pep=None):
    if pep is None: pep = sp['pep']
    # 计算肽的质量
    mass = fastmass(pep, 'M', sp['charge'], mod=sp['mod'], nmod=sp['nmod'])
    # 返回偏差值，单位为ppm
    return ((sp['mass'] - mass) / mass) * 1000000

def fastmass(pep, ion_type, charge, nmod=None, mod=None, cam=True):        
    # 计算肽段的基质量
    base = mass.fast_mass(pep, ion_type=ion_type, charge=charge)

    # 如果cam为True，则加上固定C修改的质量
    if cam: base += 57.021 * pep.count('C') / charge # fixed C modification
    
    # 返回基质量
    return base
