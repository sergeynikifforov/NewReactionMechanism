import ReactionClass
import matplotlib.pyplot as plt
h_mass =1.67e-27
temperature = [i/2. for i in range(1200, 1800)]
k1 = list()
k9_inf = list()
k9_0 = list()
k17b = list()
k_hydrogen = list()
k_hydroperoxyl_radical = list()
pressure_new = list()
temperature_new = list()
for i in range(len(temperature)):
    k1.append(ReactionClass.ReactionClass(1.04e14, 0., 15286).constant_rate_evaluate(temperature[i]))
    k9_inf.append(ReactionClass.ReactionClass(5.59e13, 0.2, 0).constant_rate_evaluate(temperature[i]))
    k9_0.append(ReactionClass.ReactionClass(2.65e19, -1.3, 0).constant_rate_evaluate(temperature[i]))
    k17b.append(ReactionClass.ReactionClass(7.15e4, 2.54, 21404).constant_rate_evaluate(temperature[i]))
    k_hydrogen.append(ReactionClass.ReactionClassWall(1e-3, 1.*h_mass).constant_rate_evaluate_wall(temperature[i]))
    k_hydroperoxyl_radical.append(ReactionClass.ReactionClassWall(1e-3, 33.*h_mass).constant_rate_evaluate_wall(temperature[i]))

for i in range(len(temperature)):
    p = ReactionClass.LimitsClass(k1[i], k17b[i], k9_inf[i], k9_0[i], k_hydrogen[i], k_hydroperoxyl_radical[i]).cubic_equation_solve(temperature[i])
    for j in range(len(p)):
        pressure_new.append(p[j])
        temperature_new.append(temperature[i])
plt.scatter(temperature_new, pressure_new, s=1)
plt.yscale('log')
plt.show()
