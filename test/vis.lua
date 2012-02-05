

local function TestTexture()

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
end


local function TestLines3d()

   local N = 100
   local points = lunum.zeros({N,3})

   for i=0,N-1 do
      local t = i/N

      local x =  0.5*(math.sin(20*t) - 0.5)
      local y =  0.5*(math.tanh(12*t) - 0.5)
      local z =  0.5*(math.cos(5*t) - 0.5)

      points[{i,0}] = x
      points[{i,1}] = y
      points[{i,2}] = z
   end
   print(points)
   visual.open_window()
   visual.draw_lines3d(points)

end


--TestTexture()
TestLines3d()

