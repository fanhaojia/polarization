from math import *
import numpy as np
import linecache
start=0;end=17;
ion_e=[10,12,6]
ne=len(ion_e);nstep=end-start+1;
n=3;m=3;mm=10;nn=[[0]*2];j=0;vec_abc={};vec_xyz={};
OUT='polar_ion.dat'
fc=open(OUT,'w')
fc=file(OUT,'a+')
ploar_ions=str(0.000)+'   '+str(0.000)+'   '+str(0.00)
fc.writelines(ploar_ions)
fc.write("\n")
for f in range(1,nstep):
    end=f;j=0;
    for i in start, end:
        fname='./CONTCAR'+str(i)+'.vasp'
        print fname
        ff=open(fname,'r'); ffline=ff.readlines();nn[0][j]=len(ffline); na=nn[0][j]+1; ff.close() #count how many lines
        print 'lines number' + '   ' +str(na)
        j=j+1;abc='abc_'+str(i);xyz='xyz_'+str(i)
        vec_abc[abc]=[[0]*m for i in range(n)]; #print abc   #latice vector
        vec_xyz[xyz]=[[0]*m for i in range(na)]; #print xyz   #atoms corordinates
        diff=[[0]*m for i in range(na)];
        vector=[[0]*m for i in range(na)];
    if nn[0][0]!=nn[0][1]:
        print "YOUR STRUCTURES HAS DIFFERENT NUMBER OF LINES"
    natom=[0]*mm; #atoms for each elements
    for i in start, end:
        fname='./CONTCAR'+str(i)+'.vasp';fp=open(fname,'r');abc='abc_'+str(i);xyz='xyz_'+str(i);
        line= fp.readline();i=1;tatom=nele=0;
        while line:
            if i in range(3,6): #latice constant
                j=0
                list=line.split(' ')
                for word in list:
                    if len(word)>0 and word != ' ' and word !='\n':#delet the '\n' at the end of line
                        if word[-1]== '\n' or word[-1]== '\t' or word[-1]== '\r' or word[-1]== '\r\n' :
                            vec_abc[abc][i-3][j]=float(word[:-1]);vec_abc[abc][i-3][j]=round(vec_abc[abc][i-3][j],16);j=j+1
                        else:
                            vec_abc[abc][i-3][j]=float(word);vec_abc[abc][i-3][j]=round(vec_abc[abc][i-3][j],16);j=j+1
            if i ==7:
                list=line.split(' ')
                j=0
                for word in list:
                    if len(word)>0 and word != ' ' and word !='\n':
                        if word[-1]== '\n' or word[-1]== '\t' or word[-1]== '\r' or word[-1]== '\r\n' :
                            natom[j]=int(word[:-1]);tatom=tatom+natom[j];j=j+1
                        else:
                            natom[j]=int(word);tatom=tatom+natom[j];j=j+1
            if i in range(9,tatom+9): #read the atoms postion(x,y,z) into matrix
                list=line.split(' ')
                j=0
                for word in list:
                    if len(word)>0 and word != ' ' and word !='\n':
                        if word[-1]== '\n' or word[-1]== '\t' or word[-1]== '\r' or word[-1]== '\r\n' :
                            vec_xyz[xyz][i][j]=float(word[:-1])
                        #xyz[i-8][j]=round(xyz[i-8][j],16)
                            j=j+1
                        else:
                            vec_xyz[xyz][i][j]=float(word[:-1])
                        #xyz[i-8][j]=round(xyz[i-8][j],16)
                            j=j+1
            line=fp.readline();i=i+1
#print vec_abc;#print vec_xyz
    fp.close()
    abc_0='abc_'+str(start);abc_1='abc_'+str(end);
    xyz_0='xyz_'+str(start);xyz_1='xyz_'+str(end);
    print 'ion_valence:' + '   ' +str(ion_e);
    print 'number_each_atoms:'+ '   ' +str(natom)
    for i in range(3):
        for j in range(m):
            diff[i][j]=float(vec_abc[abc_1][i][j]-vec_abc[abc_0][i][j])
    for i in range(9,na):
        for j in range(m):
            diff[i][j]=float(vec_xyz[xyz_1][i][j]-vec_xyz[xyz_0][i][j])
    #print diff
    natom[0]=natom[0]+9;
    for l in range(1,ne):
        natom[l]=natom[l-1]+natom[l];
    print natom
    for i in range(9,natom[0]):
        for j in range(m):
            vector[i][j]=float(diff[i][j]*ion_e[0]*1.0)
    for l in range(1,ne):
        if natom[l]<=na:
            for i in range(natom[l-1],natom[l]):
                for j in range(m):
                    vector[i][j]=float(diff[i][j]*ion_e[l]*1.0)
    #print vector
    vec_x=vec_y=vec_z=polar_x=polar_y=polar_z=0
    for i in range(9,na):
        vec_x=round(vec_x+vector[i][0],8)
        vec_y=round(vec_y+vector[i][1],8)
        vec_z=round(vec_z+vector[i][2],8)
#print vec_x, vec_y, vec_z
    polar_x=round(vec_x*vec_abc[abc_1][0][0]+vec_y*vec_abc[abc_1][1][0]+vec_z*vec_abc[abc_1][2][0],8)
    polar_y=round(vec_x*vec_abc[abc_1][0][1]+vec_y*vec_abc[abc_1][1][1]+vec_z*vec_abc[abc_1][2][1],8)
    polar_z=round(vec_x*vec_abc[abc_1][0][2]+vec_y*vec_abc[abc_1][1][2]+vec_z*vec_abc[abc_1][2][2],8)
    print polar_x, polar_y, polar_z
#new_abc=[[0]*m for i in range(n)];new_xyz=[[0]*m for i in range(na)];
    ploar_ions=str(polar_x)+'   '+str(polar_y)+'   '+str(polar_z)
    fc.writelines(ploar_ions)
    fc.write("\n")
fc.close()




















#        if i ==8:                #judge the fractional coordinates or Cartesian coordinates
#            list=line.split(' ')
#           for word in list:
#               if len(word)>0 and word != ' ':
#                   if word[0] =='d' or word[0] == 'D':
#                       head =1 #fractional coordinates
#                   elif word[0] =='c' or word[0] == 'C':
#                       head =0 #Cartesian coordinates
#                   else:
#                       print "warning the in line 7"
#       if i ==6:
#           list=line.split(' ')
#           for word in list:
#               if len(word)>0 and word != ' ':
#                   if word[-1]== '\n':
#                       element[j]=word[:-1] 
#                       j=j+1
#                   else:
#                       element[j]=word
#                       j=j+1
#           nele=j #count how many element in this POSCAR
#       if i ==7:
#           list=line.split(' ')
#           for word in list:
#               if len(word)>0 and word != ' ':
#                   if word[-1]== '\n':
#                       natom[j]=int(word[:-1]) 
#                       tatom=tatom+natom[j]
#                       j=j+1
#                   else:
#                       natom[j]=int(word)
#                       tatom=tatom+natom[j]
#                       j=j+1
