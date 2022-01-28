import prefect
# from prefect import task, Flow, Parameter
from prefect import task, Flow

@task
def add(input, y=10):
    return input + y

@task
def add2(input, y=20):
    return input + y

with Flow("additions") as f:
    # startval = Parameter("start_value")

    step1 = add(7, y=11)
    step2 = add2(step1, y=21)

    print(step2)

# f.run(start_value = 7)
state = f.run()
type(state._result.value)