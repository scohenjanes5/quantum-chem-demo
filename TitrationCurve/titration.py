# Vb(Cb+H-Kw/H) = Va(Ca*(K1*H^2 + 2*K1*K2*H + 3*K1*K2*K3)/(H^3 + K1*H^2 + K1*K2*H + K1*K2*K3)) - H + Kw/H)

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from sympy.abc import x
import argparse

def getArgs():
    parser = argparse.ArgumentParser(description='Calculate the pH of a titration')
    acid_charicteristics = parser.add_mutually_exclusive_group()
    acid_charicteristics.add_argument('-Ka', type=float, nargs='*', help='Equilibrium constants. Up to 3 allowed.', default=[1e-4,0,0])
    acid_charicteristics.add_argument('-pKa', type=float, nargs='*', help='pKa values. Up to 3 allowed.', default=[])
    parser.add_argument('-Ca', type=float, help='Concentration of the acid', default=0.1)
    parser.add_argument('-Cb', type=float, help='Concentration of the base', default=0.1)
    parser.add_argument('-Va', type=float, help='Volume of the acid stock solution', default=0.01)
    parser.add_argument('-Vw', type=float, help='Volume of water used to dilute acid.', default=0)
    parser.add_argument('-Kw', type=float, help='Ion product of water', default=1e-14)
    parser.add_argument('--export', help='Export the polynomial data to a file', action='store_true', default=False)
    parser.add_argument('--approx', help='Approximate pH for monoprotic acids', action='store_true', default=False)
    parser.add_argument('--fast', help='Use a faster method to calculate pH', action='store_true', default=False)
    args = parser.parse_args()


    if args.fast:
        args.approx = True


    num_Ka_supported_for_approx = 1
    if args.approx == True and (len(args.pKa) > num_Ka_supported_for_approx or np.count_nonzero(args.Ka) > num_Ka_supported_for_approx):
        print("approximation only works for monoprotic acids")
        quit()

    args.Ca = args.Ca*args.Va/(args.Va + args.Vw)
    args.Va = args.Va + args.Vw

    if len(args.pKa) > 0:
        Ka = [10**-pKa for pKa in args.pKa]
        args.Ka = Ka

    if len(args.Ka) < 3:
        args.Ka += [0]*(3-len(args.Ka))

    return args

def pH(args, Vb=0, name="params.txt"):

    Ka = args.Ka
    Ca = args.Ca
    Cb = args.Cb
    Kw = args.Kw
    Va = args.Va
    Na = Va*Ca
    Nb = Vb*Cb
    Vt = Va + Vb

    K12 = Ka[0]*Ka[1]
    K123 = K12*Ka[2]

    x_0 = -K123*Kw*Vt
    x_1 = -3*K123*Na + K123*Nb -K12*Kw*Vt
    x_2 = -2*K12*Na + K123*Vt -Ka[0]*Kw*Vt + Nb*K12
    x_3 = -Ka[0]*Na + K12*Vt - Kw*Vt + Ka[0]*Nb
    x_4 = Ka[0]*Vt + Nb
    x_5 = Vt

    # polynomial version of the equation
    solution = sp.solve(x_0 + x_1*x + x_2*x**2 + x_3*x**3 + x_4*x**4 + x_5*x**5, x)

    #Ignore tiny imaginary parts and select positive real roots
    solution = [sp.re(i) for i in solution if sp.im(i) < 1e-10 and sp.re(i) > 0]

    if len(solution) > 1:
        print('Multiple roots found. Cannot determine pH')
        return np.nan
    elif len(solution) == 0:
        print('No roots found. Cannot determine pH')
        return np.nan
    else:
        pH = -np.log10(float(solution[0]))

    if args.export:
        line_arr = [Ka[0], Ka[1], Ka[2], Ca, Cb, Va, Vb, x_0, x_1, x_2, x_3, x_4, x_5, pH]
        line = ','.join([str(i) for i in line_arr])
        with open(name, 'a') as f:
            f.write(line + '\n')

    #print(f'after adding {Vb}L of base to {Va}L of acid, the pH is {pH:.2f}')
    return pH

def pH_approx(args, Vb=0, use_OH=False):
    Ka = args.Ka[0]
    Ca = args.Ca
    Cb = args.Cb
    Kw = args.Kw
    Va = args.Va
    Na = Va*Ca
    Nb = Vb*Cb
    Vt = Va + Vb
    Kb = Kw/Ka

    if use_OH:
        #Use the designations in the paper
        Y5 = Kb - (Nb-Na)/Vt

        OH = 0.5 * (-Y5 + np.sqrt(Y5**2 + 4*Kb*(Nb/Vt)))
        pH = -np.log10(Kw) + np.log10(OH)
    else:
        #Use the designations in the paper
        Y3 = -Ka + Nb/Vt
        Y1 = Na-Nb

        H = 0.5 * (-Y3 + np.sqrt(Y3**2 + 4*Ka*Y1/Vt))
        pH = -np.log10(H)

    return pH

def main():
    args = getArgs()

    new_increment = 1e-4
    last_volume = 0

    Ka = args.Ka[0]
    pKa = -np.log10(Ka)

    if not args.fast:
        pH_vals = [pH(args, last_volume)]
        last_pH = pH_vals[-1]

        if args.approx:
            approx_pH = [pH_approx(args, last_volume)]
            last_pH = approx_pH[-1]
        volumes = [last_volume]
        derivative_threshold = 30
        derivative = 10

        while last_pH < 12 or derivative > derivative_threshold:
            last_volume += new_increment
            if args.approx:
                if approx_pH[-1] > pKa:
                    approx_pH.append(pH_approx(args, last_volume, use_OH=True))
                else:
                    approx_pH.append(pH_approx(args, last_volume))
                if args.fast:
                    prev_pH = last_pH
                    last_pH = approx_pH[-1]

            if not args.fast:
                pH_vals.append(pH(args, last_volume))
                prev_pH = last_pH
                last_pH = pH_vals[-1]

            volumes.append(last_volume * 1000)
            try_derivative = (last_pH - prev_pH) / new_increment
            if try_derivative != np.nan:
                derivative = try_derivative

            print(f'Volume: {last_volume:.2f}L, pH: {last_pH:.2f}, derivative: {derivative:.2f}')

    else:
        approx_pH = []

        volumes = np.linspace(0, 0.02, num=1000)
        last_pH = 0
        for vol in volumes:
            if last_pH > pKa:
                approx_pH.append(pH_approx(args, vol, use_OH=True))
            else:
                approx_pH.append(pH_approx(args, vol))
            last_pH = approx_pH[-1]
        volumes *= 1000
        approx_pH = np.array(approx_pH)

#    quit()
    if not args.export:
        if not args.fast:
            plt.plot(volumes, pH_vals, label='Exact')
        if args.approx:
            plt.plot(volumes, approx_pH, label='Approx')
        plt.title(f'{args.Ca}M acid, pKa = {float(pKa)}, Titrated by {args.Cb}M NaOH')
        plt.xlabel('Volume of base added (mL)')
        plt.ylabel('pH')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    main()
