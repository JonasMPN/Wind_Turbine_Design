import pandas as pd
save_as = "data/IEA_10MW/airfoils/airfoil_additional_information.dat"
lines_to_use = {
    "StallAngle1": 18,
    "StallAngle2": 19,
    "CnSlope": 21,
    "CritCn1": 36,
    "CritCn2": 37,
}
data = {col: list() for col in ["filename", *lines_to_use]}
for i in range(30):
    idx = f"{i}" if i>9 else f"0{i}"
    file = f"data/IEA_10MW/airfoils/FAST_format/IEA-10.0-198-RWT_AeroDyn15_Polar_{idx}.dat"
    data["filename"].append(file)
    with open(file, "r") as f:
        lines = f.readlines()
        for column, line_number in lines_to_use.items():
            value = str()
            for char in lines[line_number]:
                if char != " ":
                    value += char
                else:
                    break
            data[column].append(float(value))
pd.DataFrame(data).to_csv(save_as, index=False)