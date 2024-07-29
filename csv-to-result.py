import csv
import sys
from optimization import result

def parse_csv(filepath, ignore_not_converge = True):
    results = result.Results()
    with open(filepath, 'r') as file:
        reader = csv.reader(file, delimiter=',')

        # Read the header row to get column names
        next(reader)

        # Print each column name
        for row in reader:
            specie = row[0]
            variable = row[1]
            variation = row[2]
            converge = True if row[3] == "True" else False
            iterations = row[4]
            gradient = row[5]
            init_value = row[6]
            final_value = row[7]

            if not converge and ignore_not_converge:
                continue


            results.add_single_result(
                    result.Result(
                        specie = specie,
                        variable = variable,
                        variation = int(variation),
                        converge = converge,
                        iterations = int(iterations),
                        gradient = gradient,
                        init_value = float(init_value),
                        final_value = float(final_value),
                    ))
    return results

if len(sys.argv) != 2:
    print("Usage: python csv-to-result.py <csv_filepath>")
    sys.exit()

filepath = sys.argv[1]
results = parse_csv(filepath, False)
parsed = results.get_results_metrics()


print("\n\n--------------\n")
print("\tLatex Data - stacionaty, percent")
print("\n--------------\n")
order_conv = ["HO+HF", "R-vdW (F--H2O)", "TS", "P-vdW (HO--HF)", "F + H2O"]

for i, conv in enumerate(order_conv):
    data = parsed["convergence"][conv]
    print("\t({}, {:.2f})".format(i+1, data["conv_percent"]))

print("\n--------------\n")
print("\tLatex Data - iterations, percent")
print("\n--------------\n")
for it, data in parsed["iterations"].items():
    print("\t({}, {:.2f})".format(it, data["percent"]))
