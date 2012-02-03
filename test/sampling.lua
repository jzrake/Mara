

-- *****************************************************************************
--
-- Test of parallel sampling routines.
--
-- *****************************************************************************

local Nsamp = 10000
local verbose = false
local start = 

set_domain({0,0,0}, {1,1,1}, {16,16,16}, 8, 2)

init_prim(function(x,y,z) 
	     local r = (x^2 +y^2 + z^2)^0.5	     
	     return {1,1, 0,0,0, x/r, y/r, z/r}
	  end)

local start = os.clock()

for i=0,Nsamp do
   local x, y, z =  math.random(), math.random(), math.random()
   local P = prim_at_point{x,y,z}
   if verbose then print(lunum.array{x,y,z}, P) end
end

print(string.format("took %d samples in %f seconds",
		    Nsamp, os.clock() - start))
