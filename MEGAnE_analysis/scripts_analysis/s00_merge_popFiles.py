import pandas as pd

# === File paths ===
f1 = "../data/PPY/s01.individuals_table.txt2"   # reference file with Crop_or_Wild info
f2 = "../data/PPY/s01.individuals_table.txt"    # smaller file (missing Crop_or_Wild)
fout = "../data/PPY/s01.individuals_table.txt3"  # output merged file

# === Read both TSV files ===
df1 = pd.read_csv(f1, sep="\t")
print(df1)
df2 = pd.read_csv(f2, sep="\t")
print(df2)

# === Build a reference mapping from species name to Crop_or_Wild ===
# Example: "Pyrus_betulifolia" -> "Rootstock"
species_to_status = (
    df1[["Population", "Crop_or_Wild"]]
    .dropna()
    .drop_duplicates()
    .set_index("Population")["Crop_or_Wild"]
    .to_dict()
)
print(species_to_status)

# === Mapping from population short names in the small file to species names ===
# We'll use this to translate 'betu' -> 'Pyrus_betulifolia' etc.
pop_to_species = {
    "betu": "Pyrus_betulifolia",
    "pash": "Pyrus_pashia",
    "pyri_JP": "Pyrus_pyrifolia",
    "Sand_CN": "Pyrus_pyrifolia",
    "Sand_CN-SE": "Pyrus_pyrifolia",
    "Sand_CN-SW": "Pyrus_pyrifolia",
    "ussu": "Pyrus_ussuriensis", #wild
    "White": "Pyrus_bretschneideri" #cultivated
}
print(pop_to_species)

# === Infer Crop_or_Wild based on species name ===
# For each row in df2, map the population to its species,
# then use the reference dictionary to infer Crop/Wild status.
df2["Species_name"] = df2["Population"].map(pop_to_species)
df2["Crop_or_Wild"] = df2["Species_name"].map(species_to_status)
print(df2)

# === Replace missing statuses with empty string if inference failed ===
df2["Crop_or_Wild"] = df2["Crop_or_Wild"].fillna("")
print(df2)

# === Save the updated table ===
df2.to_csv(fout, sep="\t", index=False)

# === Summary ===
print(f"✅ Updated file written to {fout}")
print(f"Number of rows: {len(df2)}")

# === Verification preview ===
print("\nPreview of inferred Crop_or_Wild values:")
print(df2[["Population", "Species_name", "Crop_or_Wild"]].drop_duplicates().head(15))
