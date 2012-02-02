

local Nx = 100
local Ny = 240

visual.open_window()

local x = lunum.zeros{Nx,Ny}
local y = lunum.zeros{Nx,Ny}

for i,j in x:indices() do
   x[{i,j}] = 0.08 * i * math.pi
end
for i,j in y:indices() do
   y[{i,j}] = 0.04 * j * math.pi
end


for iter=0,1000 do

   local A = lunum.sin(x+math.cos(0.1*iter)) + lunum.cos(y-math.sin(0.1*iter))
   print(iter)
   collectgarbage()
   visual.draw_texture(A)

end
