

-- *****************************************************************************
--
-- Test of parallel sampling routines.
--
-- *****************************************************************************

local Nsamp = 0
local Nzone = 15
local verbose = true


local function TestSamplingNd(dims)

   if     dims == 1 then
      set_domain({0}, {1}, {Nzone}, 8, 2)

   elseif dims == 2 then
      set_domain({0,0}, {1,1}, {Nzone,Nzone}, 8, 2)

   elseif dims == 3 then
      set_domain({0,0,0}, {1,1,1}, {Nzone,Nzone,Nzone}, 8, 2)
   end

   local coords = { }
   init_prim(function(x,y,z)
		table.insert(coords, lunum.array({x,y,z}))
		if (x == 0.5 and y == 0.5) then return {1,1, 0,0,0, 0,0,0}
		else
		   return {1,1, 0,0,0, x,y,z}
		end
	     end)

   visual.draw_texture(get_prim().Bx)

   local start = os.clock()

   for i=0,Nsamp do
      local x, y, z =  math.random(), math.random(), math.random()
      local P = prim_at_point{x,y,z}
      if verbose then print(lunum.array{x,y,z}, P) end
   end
   for k,v in pairs(coords) do
      local P = prim_at_point(v)
      if verbose then print(v,P) end
   end

   print(string.format("took %d samples in %f seconds",
		       Nsamp, os.clock() - start))
end

set_fluid("rmhd")
visual.open_window()
TestSamplingNd(2)
