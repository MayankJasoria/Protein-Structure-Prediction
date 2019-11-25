from flask import Flask, render_template, request
from main import MainApp

app=Flask(__name__)

mainapp = None

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/results/", methods=["GET","POST"])
def execute():
    # run the prediction and fetch the results
    cf_prediction, gor_prediction, rost_prediction, zf_prediction, bt_prediction, dgp_prediction = mainapp.predict(request.form["sequence"])
    return render_template("results.html",
        ip_seq=request.form["sequence"],
        cf_pred=cf_prediction,
        gor_pred=gor_prediction,
        rs_pred=rost_prediction,
        bl_pred=bt_prediction,
        zf_pred=zf_prediction,
        dgp_pred=dgp_prediction
    )

if __name__ == "__main__":
    mainapp = MainApp()
    app.run(debug=True)