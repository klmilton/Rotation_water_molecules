from ase.io import read,write
from ase import Atoms
import numpy as np


def y_shift(water, ws2,config):
    wpos=water.positions
    y=[]
    
    for x in range(len(wpos)):
        y.append(wpos[x][1])
    
    wy=np.min(y)
    #print(wy)

    diff=wy-wsy
    #2.49 chosen as from geometry opt was found to be 
    #distance of water from ws2 surface
    m=2.49-diff
    water.translate([0,m,0])

    
    interface=ws2+water
    interface.write('./'+config+'_tungsten.xyz')
    
    
def shifting(water,ws2,config):
    wpos=water.positions
    x=[]
    y=[]
    z=[]
    ymin=100
    indices=0
    for ii in range(len(wpos)):
        y1=wpos[ii][1]
        #print(y1)
        x.append(wpos[ii][0])
        y.append(wpos[ii][1])
        z.append(wpos[ii][2])
        if y1 < ymin:
            ymin=y1
            indices=ii

    atom=wpos[indices]
    #numbers below are taken from the xz positions I looked at in 
    #the cell to be ontop of S and in the middle of the ring
    center=[17.290958, wpos[indices][1],16.493698]
    sulfur=[17.355930, wpos[indices][1],14.564505]
    diff_c=center-atom
    diff_s=sulfur-atom
    wc=water.copy()
    ws=water.copy()
    wc.translate([diff_c[0],diff_c[1],diff_c[2]])
    ws.translate([diff_s[0],diff_s[1],diff_s[2]])


    interface=ws2+wc
    interface.write('./'+config+'_center.xyz')
    
    interface=ws2+ws
    interface.write('./'+config+'_sulfur.xyz')  
    
    
def sort_parts(file, cl):
    #open up the xyz file required
    atoms=read(file)
    
    #set cell params here
    #atoms.set_cell([19.011854486183466,50,19.011854486183466,90,120,90])
    cell=atoms.get_cell()
    
    atoms.set_cell(cl)
    
    #can change this to oxygen if required
    water_atoms=[]
    water_atom_index=[]
    WS2_atoms=[]
    
    #loop through atoms to move into two different atoms object: WS2 and water (in my case)
    for x in atoms:
        if x.symbol == 'H' or x.symbol =='O':
            water_atoms.append(x)
            water_atom_index.append(x.index)
        else:
            WS2_atoms.append(x)
    
    
    ws2=Atoms(WS2_atoms,cell=cell,pbc=(True,False,True))
    wspos=ws2.positions
    wsy=np.max(wspos[:][1])
    water=Atoms(water_atoms,cell=cell,pbc=(True,False,True)) #pbc here set to whatever way you use - I switch my z and y axis alot. In this case my y is perpendicular to the slab

    return ws2, wsy, water

#waterpar=Atoms(water_atoms,cell=cell,pbc=(True,True,False))
#wat1=Atoms(water_atoms,cell=cell,pbc=(True,True,False))

file='initial_config.xyz'
cell=[[31.548486044253210,0.0,0.0],[-3.6432207681118972E-16,48.616204921721383E+01,0.0],[-1.5777095817208711E+01,1.1155748083267115E-13,2.7326731551312680E+01]]
    

ws2, wsy, water=sort_parts(file, cell)
ws2, wsy, waterpar=sort_parts(file, cell)
ws2, wsy, wat1=sort_parts(file, cell)
ws2, wsy, waterbothh=sort_parts(file, cell)


#the part with rotation!
water.rotate(90, 'x',center='COM',rotate_cell=False)
y_shift(water,ws2,'para')


waterpar.rotate(180, 'x',center='COM',rotate_cell=False)
y_shift(waterpar,ws2,'o_down')


y_shift(waterbothh,ws2,'both_h_down')

wat1.rotate(90, 'z',center='COM',rotate_cell=False)
y_shift(wat1,ws2,'h_down')


        
shifting(water,ws2,'para')
shifting(waterpar,ws2,'o_down')
shifting(wat1,ws2,'h_down')   
shifting(waterbothh,ws2,'both_h_down')     
        
        
        
        
        
        