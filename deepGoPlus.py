import requests
import json

class DeepGoPlus:

    def __init__(self):
        self.url = "http://deepgoplus.bio2vec.net/deepgo/api/create"
        self.mediatype={"Content-Type":"application/json"}

    def execute(self, sequence):
        payload={"data_format":"fasta",
                "data":"ALAIMYPHA",
                "threshold":0.5
                }
        response = requests.post(url=self.url, headers=self.mediatype, data=json.dumps(payload))
        #print(response.request.body)
        #print(response.request.headers)
        results = json.loads(response.text)
        results = results["predictions"][0]
        functions = results["functions"]
        scores = results["scores"]

        from itertools import chain
        c = list((zip(scores, functions)))
        c.sort()

        s = ""
        for i in range(4):
            s = s + (c[i][1] + ', ')
        s = s + c[4][1]

        return s