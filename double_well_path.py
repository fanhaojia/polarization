import numpy as np
import os
import yaml
from numpy.linalg import inv
import linecache
from pymatgen.core.periodic_table import Element



def normalize_atom_position(a):
    if a < 0.0:
        a=a+1
    elif a > 1.0:
        a=a-1
    else:
        a = a
    return a

def split_line_str(line):
    out=[]
    wlist=line.split(' ')
    for word in wlist:
        if len(word)>0 and word != ' ' :
            if word[-1]== '\n':
                if len(word[:-1])>0 and word[:-1]!='':
                    out.append(word[:-1])
            else:
                out.append(word)
    return out 

def split_line_int(line):
    out=[]
    wlist=line.split(' ')
    for word in wlist:
        if len(word)>0 and word != ' ' :
            if word[-1]== '\n':
                if len(word[:-1])>0 and word[:-1]!='':
                    out.append(int(word[:-1]))
            else:
                out.append(int(word))
    return out 

def split_line_float(line):
    out=[]
    wlist=line.split(' ')
    for word in wlist:
        if len(word)>0 and word != ' ' :
            if word[-1]== '\n':
                if len(word[:-1])>0 and word[:-1]!='':
                    out.append(float(word[:-1]))
            else:
                out.append(float(word))
    return out 

def fractional_to_cartesian(xyz, abc):
    nxyz=np.dot(xyz, abc)
    return nxyz

def cartesian_to_fractional(xyz, abc):
    nxyz=np.dot(xyz, inv(abc))
    return nxyz

def distance_two_atoms(a,b):
    dis=np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)
    return dis


def read_poscar(fname='POSCAR', cartesian=False):
    scaling=1.0
    xline = linecache.getline(fname, 6)  ##count how many element in this POSCAR
    elements=split_line_str(xline)
    #print (elements)
    xline = linecache.getline(fname, 7)  ##count how many atoms for each element in this POSCAR
    natom=split_line_int(xline)
    #print (natom) 
    xline = linecache.getline(fname, 8)  ##POSCAR is in Cartesian or fractional
    poscar_type=xline[0]
    len_block=sum(natom)+9
    abc=[];xyz=[]
    for i in range (len_block):
        if i in range(3,6):
            xline = linecache.getline(fname, i)
            abc.append(split_line_float(xline))
        elif i > 8:
            xline = linecache.getline(fname, i)
            a=split_line_float(xline)
            if poscar_type=='D':
                for h in range(3):
                    a[h]=normalize_atom_position(a[h])
                    h=h+1
                if cartesian==True:
                    b=fractional_to_cartesian(a, abc)
                else:
                    b=a
            else:
                if cartesian==True:
                    b=a
                else:
                    b=cartesian_to_fractional(a, abc)
            xyz.append(b)
        i=i+1
    
    return natom, elements, np.array(abc), np.array(xyz)

def write_poscar(natom, elements, abc, xyz, contcar='POSCAR_new.vasp', mode='Fractional'):
    tatom=sum(natom)
    str_xyz=[]
    for i in range(tatom):
        x="%.16f" % (float(xyz[i][0]))
        y="%.16f" % (float(xyz[i][1]))
        z="%.16f" % (float(xyz[i][2]))
        sstr='    '+str(x)+'    '+str(y)+'    '+str(z)
        str_xyz.append(sstr)
    a1="%.16f" % (float(abc[0][0]));a2="%.16f" % (float(abc[0][1]));a3="%.16f" % (float(abc[0][2]))
    b1="%.16f" % (float(abc[1][0]));b2="%.16f" % (float(abc[1][1]));b3="%.16f" % (float(abc[1][2]))
    c1="%.16f" % (float(abc[2][0]));c2="%.16f" % (float(abc[2][1]));c3="%.16f" % (float(abc[2][2]))
    lattice_a='    '+str(a1)+'    '+str(a2)+'    '+str(a3)
    lattice_b='    '+str(b1)+'    '+str(b2)+'    '+str(b3)
    lattice_c='    '+str(c1)+'    '+str(c2)+'    '+str(c3)
    fc=open(contcar,'w')
    fc=open(contcar,'a+')
    fc.writelines('Have added Displacements to this new POSCAR')
    fc.write("\n")
    fc.writelines(' 1.0')
    fc.write("\n")
    fc.writelines(lattice_a)
    fc.write("\n")
    fc.writelines(lattice_b)
    fc.write("\n")
    fc.writelines(lattice_c)
    fc.write("\n")
    str_na='    '; str_ele='    '
    for i in range(len(natom)):
        str_na=str_na+str(natom[i])+'    '
        str_ele=str_ele+str(elements[i])+'    '
    fc.writelines(str_ele)
    fc.write("\n")
    fc.writelines(str_na)
    fc.write("\n")
    if mode == 'Fractional':
        fc.writelines('Direct')
        fc.write("\n")
    elif mode == 'Cartesian':
        fc.writelines('Cartesian')
        fc.write("\n")     
    for i in range(tatom):
        fc.writelines(str_xyz[i])
        fc.write("\n")
    fc.close()


def get_Q(ini='POSCAR_ini', fin='POSCAR_fin'):
    natom1, elements1, abc1, xyz1=read_poscar(fname=ini, cartesian=False)
    natom2, elements2, abc2, xyz2=read_poscar(fname=fin, cartesian=False)
    disp=xyz2-xyz1
    natoms=len(xyz1)
    #get displacements in Cartesian
    disp_cart=np.dot(disp, abc1)

    ele=[]
    for i in range(len(natom1)):
        ele_i=[elements1[i] for j in range(natom1[i])]
        ele.extend(ele_i)
    mass=[Element(s).atomic_mass for s in ele]
    
    Delta_Q = 0.0
    Delta_R = 0.0
    for i in range(natoms):
        for j in range(3):
            Delta_R += (disp_cart[i][j])**2
            Delta_Q += (disp_cart[i][j])**2 * mass[i]
    R = Delta_R**0.5
    Q = Delta_Q**0.5
    print ("Amplitude Q (in angstrom*sqrt(amu)): ", Q)
    return Q

natom1, elements1, abc1, xyz1=read_poscar(fname='ini', cartesian=False)
natom2, elements2, abc2, xyz2=read_poscar(fname='fin', cartesian=False)

natoms=len(xyz1)
disp=xyz2-xyz1
#check the max displacement
dis_max=0.0
i_max=0
for i in range(natoms):
    for j in range(3):
        if disp[i][j]>dis_max:
            dis_max=disp[i][j]
            i_max=i
if dis_max>0.05:
    print ('WARNING: you have a large displacement, be carefully when checking the atom index of POSCAR_ini and POSCAR_fin, they need be kept the same')
print ("The largest displacement is (fractional): ", dis_max)
print ("It is by atom ", i_max+1, ", which is (fractional)", disp[i_max])
print ("It is by atom ", i_max+1, ", which is (Cartesian)", np.dot(disp[i_max],abc1))
#we have get the displacements matrix

nstep=10
for i in range(nstep+1):
    xyz_new=[]
    for j in range(natoms):
        d=[]
        for k in range(3):
            d.append(xyz1[j][k]+(i/nstep)*disp[j][k])
        xyz_new.append(d)
    cont='POSCAR_'+str(i)
    write_poscar(natom1, elements1, abc1, xyz_new, contcar=cont, mode='Fractional')
    
nstep2=5
for i in range(1, nstep2+1):
    xyz_new=[]
    for j in range(natoms):
        d=[]
        for k in range(3):
            d.append(xyz1[j][k]+(i+nstep2)/nstep2*disp[j][k])
        xyz_new.append(d)
    cont='POSCAR_'+str(i+nstep)
    write_poscar(natom1, elements1, abc1, xyz_new, contcar=cont, mode='Fractional')

dQ=[]
fc=open('Q.dat','w')
for i in range (0, nstep+nstep2):
    fname='POSCAR_'+str(i)
    Q_i=get_Q(ini='POSCAR_0', fin=fname)
    dQ.append(Q_i)
    fc.writelines(str(i)+'    '+str(Q_i))
    fc.write("\n")
fc.close()

