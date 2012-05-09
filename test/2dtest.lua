
local tests = require 'tests'

function RunSimulation(CFL, EndTime)

   visual.open_window()

   local dt = 0.0
   local Iteration = 0
   local CurrentTime = 0.0
   local piter = true
   local start = os.clock()

   print_mara()
   init_prim(Explosion)--KelvinHelmoltz)

--   local problem = tests.DensityWave
--   problem.velocity = { 1.0, 1.0, 0 }
--   problem.mode = { 1.0, 1.0, 0 }
--   init_prim(tests.DensityWave:get_pinit(0.0))

   while CurrentTime < EndTime do

      local prim = get_prim()
      local draw_array = prim.rho[':,:']
      visual.draw_texture(draw_array)

--      local ppmname = string.format("images/output-%04d.ppm", Iteration)
--      write_ppm(ppmname, prim["rho"])
--      os.execute("convert -rotate 90 " .. ppmname .. " " .. string.gsub(ppmname, ".ppm", ".png"))
--      os.execute("rm " .. ppmname)

      local kzps = advance(dt)

      if piter then
         print(string.format("%05d: t=%5.4f, dt=%5.4e %3.2fkz/s",
                             Iteration, CurrentTime, dt, kzps))
      end

      dt = get_timestep(CFL)
      CurrentTime = CurrentTime + dt
      Iteration = Iteration + 1

      collectgarbage()
   end

   print(string.format("run took %3.2f seconds", os.clock() - start))
   return get_prim()
end


function Explosion(x,y,z)
   local r2 = x*x + y*y + z*z
   if r2 < 0.005 then
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


--set_domain({-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5}, {16,16,16}, 5, 7)
set_domain({-0.5, -0.5}, {0.5, 0.5}, {16, 16}, 5, 3)
set_fluid("euler")
set_eos("gamma-law", 1.4)
--set_boundary("reflect2d", 2, 3)
set_boundary("periodic")
set_riemann("hllc")
--set_advance("single")

set_advance("rk4")
--set_godunov("plm-split")
set_godunov("weno-split")
RunSimulation(0.4, 6.0)

