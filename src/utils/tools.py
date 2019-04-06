

#
# Selects rows with values higher than "value" in column "column".
#
def select_data(source, column, value):
    selected_data = source[source[column] >= value]
    print(f"* Selected {len(selected_data.index)} rows.")

    return selected_data
