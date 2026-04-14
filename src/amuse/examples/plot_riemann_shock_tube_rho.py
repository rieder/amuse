from matplotlib import pyplot as plt
from pathlib import Path

import csv


data_dir = Path(__file__).parent


def read_csv(filename):
    ifile  = open(filename, "r")
    reader = csv.reader(ifile)

    x = []
    rho = []
    rownum = 0
    for row in reader:
        # Save header row.
        if rownum == 0:
            header = row
        elif rownum <=2 :
            units_str = row
        else:
            colnum = 0
            x.append(float(row[0]))
            rho.append(float(row[1]))
            for col in row:
                print('%-8s: %s' % (header[colnum], col))
                colnum += 1
        rownum += 1

    ifile.close()
    return x, rho

def plot_riemann_shock_tube_rho():

    figure = plt.figure()
    plt.xlabel("[length]")
    plt.ylabel("[mass/length$^3$]")
        
    x, rho = read_csv(data_dir / "riemann_shock_tube_problem_exact.csv")
    plt.plot(x,rho)
    # x, rho = read_csv(data_dir / "riemann_shock_tube_rho_fiN7.csv")
    # plt.scatter(x, rho, s=100, marker="o", lw=0)
    x, rho = read_csv(data_dir / "riemann_shock_tube_problem_athenaN2.csv")
    plt.scatter(x, rho, s=100, marker="s", lw=0)

    plt.xlim(0.2,0.8)

#        pyplot.savefig("riemann_shock_tube_rho_"+model.name_of_the_code+".png")
    plt.savefig("riemann_shock_tube_rho")
    plt.show()


if __name__ == "__main__":
    plot_riemann_shock_tube_rho()
    
