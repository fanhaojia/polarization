from math import *
import numpy as np
import linecache
start=0;end=3;
nstep=end-start+1
n=3;m=3;mm=10;nn=[[0]*2];j=0;vec_abc={};vec_xyz={};
for i in start, end:
    fname='./POSCAR'+str(i)+'.vasp'
    print fname
    ff=open(fname,'r'); ffline=ff.readlines();nn[0][j]=len(ffline); na=nn[0][j]+1; ff.close() #count how many lines
    j=j+1;abc='abc_'+str(i);xyz='xyz_'+str(i)
    vec_abc[abc]=[[0]*m for i in range(n)]; #print abc   #latice vector
    vec_xyz[xyz]=[[0]*m for i in range(na)]; #print xyz   #atoms corordinates
if nn[0][0]!=nn[0][1]:
    print "YOUR STRUCTURES HAS DIFFERENT NUMBER OF LINES"
natom=[0]*mm; #atoms for each elements
for i in start, end:
    fname='./POSCAR'+str(i)+'.vasp';fp=open(fname,'r');abc='abc_'+str(i);xyz='xyz_'+str(i);
    line= fp.readline();i=1;tatom=nele=0;
    while line:
    while line:
        if i ==2:
            j=0
            list=line.split(' ')
            print list
            for word in list:
                if len(word)>0 and word != ' ' and word !='\n':#delet the '\n' at the end of line
                    if word[-1]== '\n' or word[-1]== '\t' or word[-1]== '\r' or word[-1]== '\r\n' :
                        chang=float(word[:-1])*1.0;j=j+1
                    else:
                        chang=float(word)*1.0;j=j+1
        if i in range(3,6): #latice constant
            j=0
            list=line.split(' ')
            for word in list:
                if len(word)>0 and word != ' ' and word !='\n':#delet the '\n' at the end of line
                    if word[-1]== '\n' or word[-1]== '\t' or word[-1]== '\r' or word[-1]== '\r\n' :
                        vec_abc[abc][i-3][j]=float(word[:-1])*chang;vec_abc[abc][i-3][j]=round(vec_abc[abc][i-3][j],16);j=j+1
                    else:
                        vec_abc[abc][i-3][j]=float(word)*chang;vec_abc[abc][i-3][j]=round(vec_abc[abc][i-3][j],16);j=j+1
        if i ==7:
            list=line.split(' ')
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
#print vec_abc;print vec_xyz
fp.close()
diff=[[0]*m for i in range(na)]
abc_0='abc_'+str(start);abc_1='abc_'+str(end);
xyz_0='xyz_'+str(start);xyz_1='xyz_'+str(end);
for i in range(3):
    for j in range(m):
        diff[i][j]=vec_abc[abc_1][i][j]-vec_abc[abc_0][i][j]
for i in range(9,na):
    for j in range(m):
        diff[i][j]=vec_xyz[xyz_1][i][j]-vec_xyz[xyz_0][i][j]
#print diff
new_abc=[[0]*m for i in range(n)];new_xyz=[[0]*m for i in range(na)];
lie=[];strnxyz=sstrnxyz=''
for line in open(fname):
    line=line.replace('\n','').split(",")
    lie.append(line)
for i in range(nstep):
    CONTCAR='CONTCAR'+str(i)+'.vasp'
    fc=open(CONTCAR,'w')
    fc=file(CONTCAR,'a+')
    fc.writelines(lie[0])
    fc.write("\n")
    fc.writelines(lie[1])
    fc.write("\n")
    for j in range(3):
        for k in range(3):
            new_abc[j][k]="%.16f" %(vec_abc[abc_1][j][k]+i*diff[j][k]/(nstep-1))
    stra='    '+str(new_abc[0][0])+'   '+str(new_abc[0][1])+'   '+str(new_abc[0][2])
    strb='    '+str(new_abc[1][0])+'   '+str(new_abc[1][1])+'   '+str(new_abc[1][2])
    strc='    '+str(new_abc[2][0])+'   '+str(new_abc[2][1])+'   '+str(new_abc[2][2])
    fc.writelines(stra)
    fc.write("\n")
    fc.writelines(strb)
    fc.write("\n")
    fc.writelines(strc)
    fc.write("\n")
    fc.writelines(lie[5])
    fc.write("\n")
    fc.writelines(lie[6])
    fc.write("\n")
    fc.writelines('Direct')
    fc.write("\n")
    for j in range(9,na):
        for k in range(3):
            new_xyz[j][k]="%.16f" %(vec_xyz[xyz_1][j][k]+i*diff[j][k]/(nstep-1)) 
    for i in range(9,na):
        sstrnxyz='  '+str(new_xyz[i][0])+'  '+str(new_xyz[i][1])+'  '+str(new_xyz[i][2])
        fc.writelines(sstrnxyz)
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
