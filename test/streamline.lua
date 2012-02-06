

local Nzone = 32

set_fluid("euler")
set_domain({0,0,0}, {1,1,1}, {Nzone,Nzone,Nzone}, 5, 2)
init_prim(function(x,y,z)
	     return {1,1, 0.1,0.1,0.1, x,y,z}
	  end)

set_boundary("outflow")
boundary.ApplyBoundaries()

local S = streamline({0.5,0.5,0.5}, 0.4, 1e-2)
print(S)
print(S:shape('array'))

