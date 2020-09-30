import pyrefinebio

r1 = pyrefinebio.ComputationalResult.search()

for i in r1:
    print(i.id)

r = pyrefinebio.Experiment.get("SRP111833")
print("done")
