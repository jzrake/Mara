

local Nx = 100
local Ny = 240

local x = 0.08 * lunum.range(Nx) * math.pi
local y = 0.04 * lunum.range(Ny) * math.pi

local A = lunum.zeros{Nx,Ny}

for i,j in A:indices() do
   A[{i,j}] = math.sin(x[i]) + math.cos(y[j])
end

visual.open_window()
visual.draw_texture(A)
