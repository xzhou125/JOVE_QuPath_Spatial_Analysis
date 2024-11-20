// Script for remapping cluster and phenotype classifications back to mIF image for quality control and visualization in QuPath
// IMPORTANT: Ensure all labels under the "Image" heading of the csv file match the image name in QuPath

// Path to CSV file
def path = "..."
// Ensure that datapath uses backslash (\) instead of forwardslash (/)

// Color separator
def delim = ","

// Get a map from cell ID -> cell
def cells = getCellObjects()
def cellsById = cells.groupBy(c -> c.getID().toString())

// Read lines
def lines = new File(path).readLines()
def header = lines.pop().split(delim)

// Handle each line
for (def line in lines) {
    def map = lineToMap(header, line.split(delim))
    def id = map['Object ID']
    def cell = cellsById[id]
    if (cell == null) {
        println "WARN: No cell found for $id"
        continue
    } else if (cell.size() != 1) {
        println "WARN: ${cell.size()} cells for $id - will skip"
        continue
    }
    // Pull labels from appropriate headings in phenotyped csv file
    cell[0].name = [map['cluster']]
    cell[0].classifications = [map['phenotype']]
}

// Helper function to create a map from column headings -> values
Map lineToMap(String[] header, String[] content) {
    def map = [:]
    if (header.size() != content.size()) {
        throw new IllegalArgumentException("Header length doesn't match the content length!")
    }
    for (int i = 0; i < header.size(); i++)
        map[header[i]] = content[i]
    return map
}