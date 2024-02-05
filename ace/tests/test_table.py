# import pandas as pd
#     from retrieve import get_fastas

#     filename = "~/Desktop/projects/ace/data/processed/tomato.csv"
#     data = pd.read_csv(filename)
#     accs = data['Protein'].apply(lambda x: x.split('|')[1]).unique()
#     fastas = get_fastas(accs)
#     accs = generate_table(fastas)
#     print(accs)