from pyrefinebio.http import get


def generator_from_pagination(response, T):
    more_results = True

    while more_results:
        for result in response["results"]:
            yield T(**result)

        if response["next"] == None:
            more_results = False
        else:
            response = get(response["next"])
