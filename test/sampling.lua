

set_domain({0,0,0}, {1,1,1}, {16,16,16}, 8, 2)
init_prim(function(x,y,z) return lunum.zeros(8)+1 end)

local P = sampling.sample({0.1, 0.1, 0.1})
print(P)

for i=0,1000 do

end
