


function RunSimulation(CFL, EndTime)

   local dt = 0.0
   local Iteration = 0
   local CurrentTime = 0.0
   local piter = true
   local start = os.clock()

   print_mara()
   init_prim(KelvinHelmoltz)

   while CurrentTime < EndTime do

      local prim = get_prim()
      local ppmname = string.format("images/output-%04d.ppm", Iteration)
      write_ppm(ppmname, prim["rho"])
      os.execute("convert -rotate 90 " .. ppmname .. " " .. string.gsub(ppmname, ".ppm", ".png"))
      os.execute("rm " .. ppmname)

      local kzps = advance(dt)

      if piter then
         print(string.format("%05d: t=%5.4f, dt=%5.4e %3.2fkz/s",
                             Iteration, CurrentTime, dt, kzps))
      end

      dt = get_timestep(CFL)
      CurrentTime = CurrentTime + dt
      Iteration = Iteration + 1

   end
   print(string.format("run took %3.2f seconds", os.clock() - start))
   return get_prim()
end


function Explosion(x,y,z)
   local r2 = x*x + y*y
   if r2 < 0.0005 then
      return { 1.000, 1.0, 0, 0, 0 }
   else
      return { 0.125, 0.1, 0, 0, 0 }
   end
end


function KelvinHelmoltz(x,y,z)

   if math.abs(y) > 0.25 then
      rho = 1.0
      vx = -0.5
   else
      rho = 2.0
      vx =  0.5
   end

   vx = 0.02*(math.random() - 0.5) + vx
   vy = 0.02*(math.random() - 0.5)

   return { rho, 2.5, vx, vy, 0.0 }
end


function ExplosionRmhd(x,y,z)
   local r2 = x*x + y*y
   if r2 < 0.025 then
      return { 1.000, 1.0, 0, 0, 0, 4, 0, 0 }
   else
      return { 0.125, 0.1, 0, 0, 0, 4, 0, 0 }
   end
end


--set_visual()
set_domain({-0.5, -0.5}, {0.5, 0.5}, {256, 256}, 5, 3)
set_fluid("euler")
set_eos("gamma-law", 1.4)
--set_boundary("reflect2d", 2, 3)
set_boundary("outflow", 2, 3)
set_riemann("hllc")
set_advance("single")

--set_godunov("weno-split")
set_godunov("plm-muscl", 2.0, 0)
RunSimulation(0.6, 1.0)

