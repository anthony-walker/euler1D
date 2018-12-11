#Testing script
import fluid_domain as fd


nN = fd.node((1,1,1,1),1,0,0,0)
print(nN.getValues())
nN = nN+(5,10,15,20)
print(nN.getValues())
