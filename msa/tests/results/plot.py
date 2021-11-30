import matplotlib.pyplot as plt


def plot_time(fname, title, xlabel, ylabel):
    with open(fname, "r") as f:
        lines = f.readlines()
    for i in range(0, len(lines), 4):
        rt = [float(t) for t in lines[i+1].strip("\n []").split(",")]
        x = list(range(len(rt)))
        plt.plot(x, rt, label=lines[i].strip("\n "))
    plt.legend(loc="upper left")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(title)


def plot_score(fname, title, xlabel, ylabel):
    with open(fname, "r") as f:
        lines = f.readlines()
    for i in range(0, len(lines), 4):
        s = [int(t) for t in lines[i+2].strip("\n []").split(",")]
        x = list(range(len(s)))
        plt.plot(x, s, label=lines[i].strip("\n "))
    plt.legend(loc="upper left")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(title)


if __name__ == "__main__":
    # plot_time("result_m_10.txt", "running time vs number of sequences", "n: number of sequences", "running time")
    plot_time("result_n_10.txt", "running time vs length of sequences", "m: length of sequences", "running time")
    # plot_score("result_m_10.txt", "score vs number of sequences", "n: number of sequences", "sp score")
    # plot_score("result_n_10.txt", "score vs length of sequences", "m: length of sequences", "sp score")

