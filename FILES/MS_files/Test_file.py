from lxml import etree

def modify_xml_lxml(file_path, output_path):
    # Parse the XML file using lxml
    parser = etree.XMLParser(encoding="ISO-8859-1", recover=True)#to handle malformed or broken xml files
    
    try:
        tree = etree.parse(file_path, parser)
        root = tree.getroot()
        
        # Check if root is valid
        if root is None:
            print("Error: XML file is malformed or has no root element.")
            return
        
        # Modify the text
        for elem in root.iter():
            if elem.attrib.get("value") == "Agilent instrument model":
                elem.attrib["value"] = "Modified Instrument Model"

        # Save modified file
        tree.write(output_path, encoding="ISO-8859-1", xml_declaration=True)
        print(f"File modified and saved to: {output_path}")

    except etree.XMLSyntaxError as e:
        print(f"XML parsing failed: {e}")
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")

# Example Usage:
#file_path=r"C:\Users\maipa\repos\oscilations_TFG_repository\MS_files"
#modify_xml_lxml(file_path+"\TEST_file.mzXML", file_path+"\TEST_file_modified_lxml.mzXML")



from pyopenms import MzXMLFile, MSExperiment

# Load and modify XML with pyOpenMS
def modify_xml_pyopenms(file_path, output_path):
    exp = MSExperiment()
    try:
        MzXMLFile().load(file_path, exp)
        print("mzXML file loaded successfully with pyOpenMS.")
    except Exception as e:
        print("Error loading mzXML with pyOpenMS:", e)
        return

    # Opens and modifies XML content as plain text
    with open(file_path, 'r', encoding="ISO-8859-1") as file:
        xml_content = file.read()

    # Replaces the "name=Conversion to mzML" with "Conversion to Modified Format"
    xml_content = xml_content.replace('name="Conversion to mzML"', 'name="Conversion to Modified Format"')

    # Save the modified XML
    with open(output_path, 'w', encoding="ISO-8859-1") as file:
        file.write(xml_content)

    print("File modified using pyOpenMS and saved to:", output_path)

# Example Usage:
file_path=r"C:\Users\maipa\repos\oscilations_TFG_repository\MS_files"
modify_xml_pyopenms(file_path+"\TEST_file.mzXML", file_path+"\TEST_file_modified_pyopenms.mzXML")
