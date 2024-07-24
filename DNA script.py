import pandas as pd
import os
import xml.etree.ElementTree as ET
import numpy as np 
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
import threading
import sys



def get_file_paths(folder_path, diffrent_folder=False):
    if not diffrent_folder:
        out_folder_path = os.path.dirname(folder_path)
    else:
        out_folder_path = folder_path
    data_file_path = os.path.join(out_folder_path,'sequencing_summary.csv')
    matches_file_path = os.path.join(out_folder_path, 'DNA_matches.csv')
    settings_file_path = os.path.join(out_folder_path, 'settings.csv')
    return data_file_path,matches_file_path,settings_file_path


# Function to extract header values
def extract_header_info(file_content):
    case_id = ""
    reading_by = ""
    reading_datetime = ""

    for line in file_content:
        if "Project:" in line:
            case_id = line.split("\\")[-2]  # Assumes case ID is in a specific position in the path
        if "Software Package:" in line:
            reading_by = line.split(":")[1].strip()
        if "Date/Time:" in line:
            reading_datetime = line.split(":")[1].strip()
    return case_id, reading_by, reading_datetime

def ensure_files_exist(folder_path, diffrent_folder=False):
    """ Ensure all necessary files exist, create them if they don't. """
    data_file_path,matches_file_path,settings_file_path=get_file_paths(folder_path,diffrent_folder=diffrent_folder)
    
    if not os.path.exists(data_file_path):
        print(f"Creating empty data file at {data_file_path}")
        pd.DataFrame(columns=['FileName','CaseID', 'SpecimenID', 'SpecimenComment', 'LocusName', 'ReadingBy', 'ReadingDateTime', 'AlleleValue']).to_csv(data_file_path, index=False)
    
    if not os.path.exists(matches_file_path):
        print(f"Creating empty matches file at {matches_file_path}")
        pd.DataFrame(columns=['LocusName', 'SpecimenID1', 'SpecimenID2', 'MatchScore', 'LatestMatchTime']).to_csv(matches_file_path, index=False)
    
    if not os.path.exists(settings_file_path):
        print(f"Creating default settings file at {settings_file_path}")
        pd.DataFrame({'Sensitivity': [0.8], 'ScannedFiles': [""]}).to_csv(settings_file_path, index=False)

def load_data(file_path):
    """ Load the consolidated data from a CSV file. """
    return pd.read_csv(file_path)

def load_existing_matches(matches_file_path):
    """ Load existing matches from a CSV file, return set of unique SpecimenIDs and DataFrame of matches. """
    matches_df = pd.read_csv(matches_file_path)
    existing_ids = set(matches_df['SpecimenID1']).union(set(matches_df['SpecimenID2']))
    return existing_ids, matches_df

def load_settings(settings_file_path):
    """ Load settings from a CSV file. """
    settings_df = pd.read_csv(settings_file_path)
    sensitivity = settings_df['Sensitivity'].iloc[0]
    scanned_files = np.unique(settings_df['ScannedFiles'].dropna().tolist())
    return sensitivity, [i for i in scanned_files], settings_df

def save_settings(settings_df, settings_file_path):
    """ Save settings to a CSV file. """
    settings_df.to_csv(settings_file_path, index=False)
def find_matches(df, sensitivity, existing_ids):
    """ Find matches between specimens based on a sensitivity threshold, excluding already matched SpecimenIDs. """
    if df.empty or 'LocusName' not in df.columns or 'ReadingDateTime' not in df.columns:
        print("No data to process or missing required columns.")
        return pd.DataFrame(columns=['LocusName', 'SpecimenID1', 'SpecimenID2', 'MatchScore', 'LatestMatchTime'])

    # Group by SpecimenID and LocusName, then aggregate alleles into a list
    grouped = df[df['SpecimenID'].apply(lambda x: x not in existing_ids)].groupby(['SpecimenID', 'LocusName']).agg({
        'AlleleValue': list,
        'CaseID': 'first',
        'SpecimenComment': 'first',
        'ReadingBy': 'first',
        'ReadingDateTime': 'first'
    }).reset_index()

    # Further group by SpecimenID to get all loci for each specimen
    specimens_grouped = grouped.groupby('SpecimenID').agg({
        'LocusName': list,
        'AlleleValue': list,
        'CaseID': 'first',
        'SpecimenComment': 'first',
        'ReadingBy': 'first',
        'ReadingDateTime': 'first'
    }).reset_index()

    matches = pd.DataFrame(columns=['SpecimenID1', 'SpecimenID2', 'MatchScore', 'LatestMatchTime'])
    matches_list = []

    specimen_ids = specimens_grouped['SpecimenID'].tolist()
    loci = specimens_grouped['LocusName'].tolist()
    alleles = specimens_grouped['AlleleValue'].tolist()
    case_ids = specimens_grouped['CaseID'].tolist()
    comments = specimens_grouped['SpecimenComment'].tolist()
    readers = specimens_grouped['ReadingBy'].tolist()
    times = specimens_grouped['ReadingDateTime'].tolist()

    n = len(specimen_ids)
    for i in range(n):
        for j in range(i + 1, n):
            common_loci = set(loci[i]).intersection(set(loci[j]))
            total_loci = len(loci[i])  # Assuming both specimens have the same number of loci
            identical_loci = 0

            for locus in common_loci:
                idx_i = loci[i].index(locus)
                idx_j = loci[j].index(locus)

                # Ensure alleles match exactly in the same order
                if alleles[i][idx_i] == alleles[j][idx_j]:
                    identical_loci += 1

            if total_loci > 0:
                match_score = identical_loci / total_loci
            else:
                match_score = 0
            if match_score >= sensitivity:
                latest_time = max(times[i], times[j])
                matches_list.append({
                    'SpecimenID1': specimen_ids[i],
                    'SpecimenID2': specimen_ids[j],
                    'MatchScore': match_score,
                    'LatestMatchTime': latest_time,
                    'CaseID1': case_ids[i],
                    'CaseID2': case_ids[j],
                    'SpecimenComment1': comments[i],
                    'SpecimenComment2': comments[j],
                    'ReadingBy1': readers[i],
                    'ReadingBy2': readers[j]
                })

    if matches_list:
        matches = pd.concat([matches, pd.DataFrame(matches_list)], ignore_index=True)

    return matches
def process_xml_file(file_path, ns):
    """ Process a single XML file and return the data as a DataFrame. """
    file_name=get_folder_file_name(file_path)
    tree = ET.parse(file_path)
    root = tree.getroot()

    # Determine format by checking the root tag or other unique identifiers
    if 'CODISImportFile' in root.tag:
        return process_codis_format(root, ns, file_name)
    elif 'DNADataTransaction' in root.tag:
        return process_niem_format(root, ns, file_name)
    else:
        print("Unknown XML format.")
        return pd.DataFrame()  # Return an empty DataFrame if format is not recognized

def process_codis_format(root, ns, file_name):
    """ Process the CODIS format XML. """
    data = []
    for specimen in root.findall('ns:SPECIMEN', ns):
        extract_common_fields(specimen, data, ns, file_name)
    return pd.DataFrame(data, columns=['CaseID', 'SpecimenID', 'SpecimenComment', 'LocusName', 'ReadingBy', 'ReadingDateTime', 'AlleleValue', 'FileName'])

def process_niem_format(root, ns, file_name):
    """ Process the NIEM format XML. """
    data = []
    for locus in root.findall('.//biom:DNALocus', ns):
        case_id = root.find('.//nc:IdentificationID', ns).text
        specimen_id = root.find('.//biom:DNASourceIdentification/nc:IdentificationID', ns).text
        specimen_comment = "Extracted from NIEM"
        locus_name = locus.find('biom:DNALocusName', ns).text
        reading_datetime = locus.find('biom:ProcessUTCDate', ns).text if locus.find('biom:ProcessUTCDate', ns) is not None else "Unknown"
        reading_by = root.find('.//biom:DeviceName', ns).text if root.find('.//biom:DeviceName', ns) is not None else "N/A"
        
        allele_values = locus.findall('.//biom:DNAAllele/biom:DNAAlleleCall1Text', ns)
        for allele in allele_values:
            allele_value = allele.text
            data.append([
                case_id, specimen_id, specimen_comment, locus_name,
                reading_by, reading_datetime, allele_value, file_name
            ])
    return pd.DataFrame(data, columns=['CaseID', 'SpecimenID', 'SpecimenComment', 'LocusName', 'ReadingBy', 'ReadingDateTime', 'AlleleValue', 'FileName'])

def extract_common_fields(specimen, data, ns, file_name):
    """ Extract fields common to both XML formats. """
    case_id = specimen.get('CASEID')
    specimen_id = specimen.find('ns:SPECIMENID', ns).text
    specimen_comment = specimen.find('ns:SPECIMENCOMMENT', ns).text if specimen.find('ns:SPECIMENCOMMENT', ns) is not None else "N/A"
    for locus in specimen.findall('ns:LOCUS', ns):
        locus_name = locus.find('ns:LOCUSNAME', ns).text
        reading_by = locus.find('ns:READINGBY', ns).text if locus.find('ns:READINGBY', ns) is not None else "N/A"
        reading_datetime = locus.find('ns:READINGDATETIME', ns).text if locus.find('ns:READINGDATETIME', ns) is not None else "N/A"
        for allele in locus.findall('ns:ALLELE', ns):
            allele_value = allele.find('ns:ALLELEVALUE', ns).text
            data.append([
                case_id, specimen_id, specimen_comment, locus_name,
                reading_by, reading_datetime, allele_value, file_name
            ])

def process_csv_file(file_path):
    df=pd.read_csv(file_path)
    file_name = get_folder_file_name(file_path)  # Extract filename from file path
    df.columns = df.columns.str.strip()
    expanded_rows = []
    for _, row in df.iterrows():
        for allele_num in range(1, 3):  # Assuming there are only two alleles max as shown
            allele_value = row[f'Allele {allele_num}'].strip()
            if allele_value:  # Ensure non-empty alleles are processed
                expanded_rows.append({
                    'CaseID': row['Sample Name'].split()[0],  # Derived from Sample Name
                    'SpecimenID': row['Sample File'],  # Using Sample File as SpecimenID
                    'SpecimenComment': file_name,  # Example comment
                    'LocusName': row['Marker'],
                    'ReadingBy': 'ABI3500',  # Static example
                    'ReadingDateTime': '0000-00-00T00:00:00',
                    'AlleleValue': allele_value,
                    'FileName': row['Sample File']  # Assuming FileName is needed
                })
    return pd.DataFrame(expanded_rows)


def scan_and_process_files(folder_path, scanned_files):
    """ Scan the folder and all subfolders for XML, TXT, and CSV files and process them into a DataFrame. """
    ns = {
        'ns': 'urn:CODISImportFile-schema',
        'biom': 'http://release.niem.gov/niem/domains/biometrics/5.1/',
        'nc': 'http://release.niem.gov/niem/niem-core/5.0/'
    }
    all_data = []
    # Walk through all directories and files in the folder path
    for root, dirs, files in os.walk(folder_path):
        for file_name in files:
            file_path = os.path.join(root, file_name)
            if file_name.endswith('.xml') and file_name not in scanned_files:
                try:
                    print(f"Found new XML file: {file_name}, processing...")
                    file_data = process_xml_file(file_path, ns)  # Existing XML processing function
                    all_data.append(file_data)
                    scanned_files.append(file_name)
                except Exception as e:
                    print(f"Error processing {file_name}: {e}")
                
            elif file_name.endswith('.txt') and file_name not in scanned_files:
                print(f"Found new TXT file: {file_name}, processing...")
                try:
                    file_data = process_txt_file(file_path)  # Your TXT file processing function
                    all_data.append(file_data)
                    scanned_files.append(file_name)
                except Exception as e:
                    print(f"Error processing {file_name}: {e}")

            elif file_name.endswith('.csv') and file_name not in scanned_files:
                print(f"Found new CSV file: {file_name}, processing...")
                try:
                    file_data = process_csv_file(file_path)  # Assume CSV files can be loaded directly
                    all_data.append(file_data)
                    scanned_files.append(file_name)
                except Exception as e:
                    print(f"Error processing {file_name}: {e}")

    if all_data:
        combined_data = pd.concat(all_data, ignore_index=True)
        print("All new files processed.")
        return combined_data, scanned_files
    else:
        print("No new files to process.")
        return pd.DataFrame(), scanned_files
    
def remove_duplicates(df, subset_columns):
    """Remove duplicates from the DataFrame based on the subset of columns."""
    original_count = len(df)
    df.drop_duplicates(subset=subset_columns, keep='last', inplace=True)
    new_count = len(df)
    print(f"Removed {original_count - new_count} duplicates; {new_count} entries remain.")
    return df


def unmelting_data(df):
    if df.empty:
        return df

    # Step 1: Clean and standardize 'LocusName' values
    df['LocusName'] = df['LocusName'].str.strip()  # Remove leading/trailing spaces in LocusName values
    locus_name_map = {
        'Amelogenin': 'AMEL',  # Mapping 'Amelogenin' to 'AMEL'
        # Add additional mappings if necessary
    }
    df['LocusName'] = df['LocusName'].replace(locus_name_map)  # Apply name changes

    # Step 2: Pivot the table to consolidate rows to a single line per specimen with loci as columns
    pivoted_data = df.pivot_table(index=['FileName', 'CaseID', 'SpecimenID', 'SpecimenComment', 'ReadingBy', 'ReadingDateTime'],
                                  columns='LocusName', values='AlleleValue', aggfunc=lambda x: list(x))

    # Step 3: Reset the index for easier data handling
    unmelted_data = pivoted_data.reset_index()

    # Step 4: Fill NaN with empty lists for uniform data handling and split lists into two sorted columns
    for locus in pivoted_data.columns:
        # Ensure all entries are lists
        filled_data = unmelted_data[locus].apply(lambda x: x if isinstance(x, list) else [])
        # Split and sort the lists into two columns, handling None values correctly
        unmelted_data[f"{locus}_1"], unmelted_data[f"{locus}_2"] = zip(
            *filled_data.apply(lambda x: sorted(x + [None, None], key=lambda v: (v is None, v))[:2])
        )

    # Step 5: Drop the original allele list columns as they are now split and sorted
    unmelted_data.drop(columns=pivoted_data.columns, inplace=True)

    # Step 6: Sort by 'ReadingDateTime' in ascending order
    unmelted_data = unmelted_data.sort_values('ReadingDateTime', ascending=True)

    return unmelted_data


def get_folder_file_name(file_path):
    """Construct a file name that includes its parent folder's name."""
    folder_name = os.path.basename(os.path.dirname(file_path))  # Get the name of the parent directory
    file_name = os.path.basename(file_path)  # Get the base file name
    return os.path.join(folder_name, file_name)  # Combine folder name and file name

# Parsing logic with automated header extraction and robust data handling
def process_txt_file(file_path):
    file_name = get_folder_file_name(file_path)  # Extract filename from file path

    with open(file_path, 'r') as file:
        file_content = file.readlines()

    case_id, reading_by, reading_datetime = extract_header_info(file_content)

    allele_data = []
    current_specimen_id = None
    data_start = False

    for line in file_content:
        if line.startswith("\tSample"):
            data_start = True
            continue
        if data_start and line.strip():
            parts = line.strip().split("\t")
            clean_parts = [part.strip() for part in parts if part.strip()]

            if clean_parts[0].isdigit():
                current_specimen_id = clean_parts[1]
                locus_name = clean_parts[2]
                allele_values = [part for part in clean_parts[3:] if part.replace('.', '', 1).isdigit() or part in ['X', 'Y']]
            else:
                locus_name = clean_parts[0]
                allele_values = [part for part in clean_parts[1:] if part.replace('.', '', 1).isdigit() or part in ['X', 'Y']]

            for allele_value in allele_values:
                allele_data.append({
                    "CaseID": case_id,
                    "SpecimenID": current_specimen_id,
                    "SpecimenComment": "",
                    "LocusName": locus_name,
                    "ReadingBy": reading_by,
                    "ReadingDateTime": reading_datetime,
                    "AlleleValue": allele_value,
                    "FileName": file_name  # Include the file name here
                })

    return pd.DataFrame(allele_data)

def main(xml_folder_path, save_to_folder_path=None, sensitivity=None):

    if save_to_folder_path is None:
        ensure_files_exist(xml_folder_path)
        data_file_path, matches_file_path, settings_file_path = get_file_paths(xml_folder_path)
        unmelted_df_file_path = os.path.join(os.path.dirname(xml_folder_path), 'final_DNA_sequencing_summary.csv')
    else:
        ensure_files_exist(save_to_folder_path, diffrent_folder=True)
        data_file_path, matches_file_path, settings_file_path = get_file_paths(save_to_folder_path, diffrent_folder=True)
        unmelted_df_file_path = os.path.join(save_to_folder_path, 'final_DNA_sequencing_summary.csv')


    if sensitivity is None:
        sensitivity, scanned_files, settings_df = load_settings(settings_file_path)
    else:
        _, scanned_files, settings_df = load_settings(settings_file_path)

    df = load_data(data_file_path)
    new_data, scanned_files = scan_and_process_files(xml_folder_path, scanned_files)  # Updated function call

    
    existing_ids, existing_matches = load_existing_matches(matches_file_path)
    # Find matches
    new_matches = find_matches(df, sensitivity, existing_ids)

    # Combine new matches with existing, sort by LatestMatchTime, and write to CSV
    if not new_matches.empty:
        full_matches = pd.concat([existing_matches, new_matches])
        full_matches['LatestMatchTime'] = pd.to_datetime(full_matches['LatestMatchTime'])
        full_matches = full_matches.sort_values('LatestMatchTime')
        full_matches.to_csv(matches_file_path, index=False)
        print("Updated matches have been saved to 'DNA_matches.csv'.")
    else:
        print("No new matches found to append.")


    if not new_data.empty:
        df = pd.concat([df, new_data], ignore_index=True)
        df=remove_duplicates(df, df.columns.tolist())
        df.to_csv(data_file_path, index=False)
        unmelted_df = unmelting_data(df)
        unmelted_df.to_csv(unmelted_df_file_path, index=False)
        print(f"Data saved to {unmelted_df_file_path}")

    # Save updated settings with scanned files

    sensitivity_series =[np.nan] * len(scanned_files)
    sensitivity_series[0]=sensitivity
  
    # Append NaNs to the DataFrame column
    settings_df={'Sensitivity':sensitivity_series,'ScannedFiles':scanned_files}
    settings_df=pd.DataFrame(settings_df)

    save_settings(settings_df, settings_file_path)

    

#================================================UI functions================================================

def run_script(input_folder_path, output_folder_path, sensitivity):
    try:
        main(input_folder_path, output_folder_path, sensitivity)
        print("Processing completed successfully!")
    except Exception as e:
        print(f"An error occurred: {e}")
        
def choose_folder(entry, other_entry):
    folder_path = filedialog.askdirectory()
    if folder_path:
        entry.delete(0, tk.END)
        entry.insert(0, folder_path)
        # Update the history file based on which entry was updated
        if entry == input_folder_path_entry:
            save_history(folder_path, other_entry.get())
        else:
            save_history(other_entry.get(), folder_path)

def start_thread():
    input_folder_path = input_folder_path_entry.get()
    output_folder_path = output_folder_path_entry.get()
    sensitivity = sensitivity_entry.get()  # Read sensitivity from entry as a string
    if os.path.exists(input_folder_path):
        save_history(input_folder_path, output_folder_path, sensitivity)  # Save the current settings to history
        thread = threading.Thread(target=run_script, args=(input_folder_path, output_folder_path, float(sensitivity)))
        thread.start()
    else:
        messagebox.showerror("Error", "The specified input folder path does not exist.")

def show_faq():
    faq_message = """
Application Functionality
1. File Processing:
    • The application will automatically scan the input folder for XML and TXT files.
    • XML files will be processed according to their format (CODIS or NIEM).
    • TXT files will be parsed to extract allele data.

2. Data Consolidation:
    • Extracted data from all files will be consolidated into a single DataFrame.
    • The data will be saved in sequencing_summary.csv in the output folder.

3. Formatting the Data:
    • The application will pivot the data to consolidate rows into a single line per specimen, with loci as columns.
    • The processed data will be saved in final_DNA_sequencing_summary.csv in the output folder.

4. Finding Matches:
    • The application will compare specimens to find matches based on the set sensitivity.
    • Existing matches will be loaded and new matches will be appended.
    • The matches will be saved in DNA_matches.csv in the output folder.

Sample Matching Methodology
1. Loading Existing Matches:
    • Existing matches from DNA_matches.csv will be loaded to ensure that previously matched specimens are not reprocessed.

2. Sensitivity Use:
    • The sensitivity value determines how stringent the matching criteria are. A higher value (closer to 1) requires more loci to match exactly.
    • Matches are found by comparing alleles for each locus. If the alleles match exactly, the locus is considered a match. Example: 24/26 alleles are identical matches, then the sensitivity score is 24/26 = 0.92.

3. Finding New Matches:
    • Specimens are grouped by LocusName, and alleles are compared to identify matches.
    • Only specimens not previously matched are considered.
    • Matches with a score above the sensitivity threshold are saved.

Notes
1. Duplicate Files:
    • If two identical copies of one file exist under different names, they will appear twice in the table and will show a full match since the table doesn't drop duplicates.

2. Resetting Result Files:
    • To reset the result files, you can delete them in the output folder or select a new output folder for a new analysis of the input folder.

Viewing Results
1. Sequencing Summary:
    • The consolidated sequencing data is saved in sequencing_summary.csv.
    • The final processed data is saved in final_DNA_sequencing_summary.csv.

2. DNA Matches:
    • Matches are saved in DNA_matches.csv.
    • The results include specimen IDs, match scores, and the latest match time.

3. Settings:
    • The settings, including the sensitivity and scanned files, are saved in settings.csv.

Technical Support
For further questions or technical support, please contact: eitanfass1996@gmail.com
"""
    messagebox.showinfo("Info", faq_message)



def load_history():
    paths = {'InputFolderPath': '', 'OutputFolderPath': '', 'Sensitivity': '0.8'}  # Default sensitivity
    try:
        with open(HISTORY_FILE_PATH, 'r', encoding='utf-8') as file:
            for line in file:
                key, value = line.strip().split('=')
                paths[key] = value.strip()
            # print("History loaded:", paths)  # Debug print
        return paths['InputFolderPath'], paths['OutputFolderPath'], paths['Sensitivity']
    except FileNotFoundError:
        print("History file not found.")  # Alert if the file is not found
        return '', '', '0.8'  # Return defaults
    except Exception as e:
        print("Failed to load history:", e)  # Print other exceptions
        return '', '', '0.8'  # Return defaults

def save_history(input_folder_path, output_folder_path, sensitivity):
    try:
        with open(HISTORY_FILE_PATH, 'w', encoding='utf-8') as file:
            file.write(f"InputFolderPath={input_folder_path}\n")
            file.write(f"OutputFolderPath={output_folder_path}\n")
            file.write(f"Sensitivity={sensitivity}\n")  # Add the sensitivity setting
        # print("History saved:", input_folder_path, output_folder_path, sensitivity)  # Debug print
    except Exception as e:
        print("Failed to save history:", e)  # Print any errors encountered

HISTORY_FILE_PATH='history.txt'

if __name__ == '__main__':
    # main(r"C:\Users\Eitan.F\OneDrive - tierraspec\Documents\איתן- אישי\DNA lab")
    root = tk.Tk()
    root.title("DNA Sequencing Data Processor")

    # Load the last used directory paths and sensitivity from history file
    input_folder_path, output_folder_path, last_used_sensitivity = load_history()

    input_folder_path_entry = tk.Entry(root, width=80)
    input_folder_path_entry.insert(0, input_folder_path)
    input_folder_path_entry.pack(pady=10)

    input_folder_select_button = tk.Button(root, text="Select Input Folder", command=lambda: choose_folder(input_folder_path_entry, output_folder_path_entry))
    input_folder_select_button.pack(pady=10)

    output_folder_path_entry = tk.Entry(root, width=80)
    output_folder_path_entry.insert(0, output_folder_path)
    output_folder_path_entry.pack(pady=10)

    output_folder_select_button = tk.Button(root, text="Select Output Folder", command=lambda: choose_folder(output_folder_path_entry, input_folder_path_entry))
    output_folder_select_button.pack(pady=10)


    sensitivity_label = tk.Label(root, text="Match Sensitivity (0-1):")
    sensitivity_label.pack()
    sensitivity_entry = tk.Entry(root, width=10)
    sensitivity_entry.insert(0, last_used_sensitivity)  # Set the last used sensitivity
    sensitivity_entry.pack(pady=5)

    start_button = tk.Button(root, text="Start Processing", command=start_thread)
    start_button.pack(pady=20)

    help_button = tk.Button(root, text="Help/Info", command=show_faq)
    help_button.pack(pady=10)

    log = scrolledtext.ScrolledText(root, state='disabled', width=70, height=10)
    log.pack(pady=20)

    def redirector(inputStr):
        log.configure(state='normal')
        log.insert(tk.END, inputStr)
        log.configure(state='disabled')
        log.yview(tk.END)

    sys.stdout.write = redirector

    root.mainloop()
