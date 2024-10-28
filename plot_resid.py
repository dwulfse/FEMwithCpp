import matplotlib.pyplot as plt

def plot_resid():
    x = []
    u = []
    u_exact = []
    residuals = []

    with open("residuals.csv", "r") as f:
        next(f)
        for line in f:
            data = line.split(",")
            x.append(float(data[0]))
            u.append(float(data[1]))
            u_exact.append(float(data[2]))
            residuals.append(float(data[3]))

    plt.plot(x, u)
    plt.plot(x, u_exact)
    plt.plot(x, residuals)
    plt.legend(["u", "u_exact", "residuals"])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("FEM approximation")
    plt.show()

def plot_errors():
    x = []
    errors = []

    with open("errors.csv", "r") as f:
        next(f)
        for line in f:
            data = line.split(",")
            x.append(float(data[0]))
            errors.append(float(data[1]))

    # log-log axis
    plt.xscale("log")
    plt.yscale("log")
    plt.plot(x, errors)
    plt.xlabel("h")
    plt.ylabel("Error")
    plt.title("Errors")
    plt.show()

plot_resid()
plot_errors()
