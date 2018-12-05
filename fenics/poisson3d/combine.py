result = open("output/result.txt", "w")

for i in range(0, 5):
    with open("output/result_"+str(i)+".txt", "r") as f:
        temp = f.read()
        result.write(temp)

result.close()
