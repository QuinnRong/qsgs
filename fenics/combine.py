import sys

result = open("output/" + sys.argv[1] + "/result.txt", "w")

num = int(sys.argv[2])
for i in range(0, num):
    with open("output/" + sys.argv[1] + "/result_"+str(i)+".txt", "r") as f:
        temp = f.read()
        result.write(temp)

result.close()
