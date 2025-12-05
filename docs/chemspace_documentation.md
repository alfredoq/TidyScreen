# ChemSpace Class Documentation

## Overview

The `ChemSpace` class is designed to manage chemical compound data within a TidyScreen project. It inherits from `ActivateProject` to access project information and provides functionality to read CSV files containing SMILES strings, compound names, and flags, storing them in a dedicated SQLite database within the project directory.

## Class Definition

```python
class ChemSpace:
    """
    ChemSpace class for managing chemical compound data within a project.
    Uses an ActivateProject object to access project information and database functionality.
    """
```

## Features

### 1. **Automatic Database Setup**
- Creates a `chemspace.db` SQLite database in the project's `data/` directory
- Automatically initializes the database schema with proper indexing
- Handles database creation and table setup transparently

### 2. **CSV File Processing** 
- Reads CSV files containing chemical compound data
- Configurable column mapping for SMILES, names, and flags
- Handles duplicate detection and skipping
- Provides detailed feedback on import results

### 3. **Data Management**
- Store compounds with SMILES strings, names, flags, and metadata
- Search compounds by name, SMILES, or flag
- Filter and retrieve compounds with various criteria
- Export data back to CSV format

### 4. **Database Operations**
- Full CRUD operations on compound data
- Optimized queries with proper indexing
- Transaction safety and error handling
- Database statistics and reporting

## Database Schema

The ChemSpace database uses the following table structure:

```sql
CREATE TABLE compounds (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    smiles TEXT NOT NULL,
    flag TEXT,
    created_date TEXT NOT NULL,
    source_file TEXT,
    UNIQUE(name, smiles)
);
```

**Indexes:**
- `idx_compounds_name` on `name` column
- `idx_compounds_smiles` on `smiles` column  
- `idx_compounds_flag` on `flag` column

## Usage Examples

### Basic Usage

```python
import tidyscreen
from chemspace.chemspace import ChemSpace

# First, activate an existing project
project_obj = tidyscreen.ActivateProject("my_chemistry_project")

# Initialize ChemSpace with the project object
chemspace = ChemSpace(project_obj)

# Load compounds from CSV file
result = chemspace.load_csv_file(
    csv_file_path="compounds.csv",
    smiles_column="smiles",
    name_column="compound_name", 
    flag_column="status"
)

if result['success']:
    print(f"Loaded {result['compounds_added']} compounds")
```

### Advanced Usage

```python
# Search for specific compounds
aspirin_compounds = chemspace.search_compounds("aspirin", "name")

# Get compounds with specific flag
approved_drugs = chemspace.get_compounds(flag_filter="approved", limit=10)

# Export filtered data
chemspace.export_to_csv("approved_compounds.csv", flag_filter="approved")

# Get database statistics
stats = chemspace.get_database_stats()
print(f"Total compounds: {stats['total_compounds']}")

# Print formatted summary
chemspace.print_database_summary()
```

## Method Reference

### Constructor

#### `__init__(project_obj: ActivateProject)`
Initializes ChemSpace with an ActivateProject object.

**Parameters:**
- `project_obj`: An instantiated ActivateProject object

**Raises:**
- `TypeError`: If project_obj is not an ActivateProject instance
- `ValueError`: If project doesn't exist in the database

### Data Loading

#### `load_csv_file(csv_file_path, smiles_column='smiles', name_column='name', flag_column='flag', skip_duplicates=True)`
Loads compounds from a CSV file into the database.

**Parameters:**
- `csv_file_path`: Path to the CSV file
- `smiles_column`: Column name containing SMILES strings
- `name_column`: Column name containing compound names
- `flag_column`: Column name containing flags
- `skip_duplicates`: Whether to skip duplicate compounds

**Returns:**
- Dictionary with keys: `success`, `message`, `compounds_added`, `duplicates_skipped`, `errors`

### Data Retrieval

#### `get_compounds(limit=None, flag_filter=None, name_pattern=None)`
Retrieves compounds with optional filtering.

**Parameters:**
- `limit`: Maximum number of compounds to return
- `flag_filter`: Filter by specific flag value
- `name_pattern`: Filter by name pattern (SQL LIKE)

**Returns:**
- List of compound dictionaries

#### `search_compounds(search_term, search_in='name')`
Searches for compounds by name, SMILES, or flag.

**Parameters:**
- `search_term`: Term to search for
- `search_in`: Field to search in ('name', 'smiles', 'flag')

**Returns:**
- List of matching compound dictionaries

#### `get_compound_count(flag_filter=None)`
Returns the total number of compounds, optionally filtered by flag.

### Data Export

#### `export_to_csv(output_path, flag_filter=None)`
Exports compounds to a CSV file.

**Parameters:**
- `output_path`: Path for the output CSV file
- `flag_filter`: Optional filter by flag value

**Returns:**
- Boolean indicating success

### Database Management

#### `get_database_stats()`
Returns comprehensive database statistics.

**Returns:**
- Dictionary with total compounds, flag counts, file counts, date ranges

#### `print_database_summary()`
Prints a formatted summary of the database contents.

#### `clear_database(confirm=True)`
Clears all compounds from the database.

**Parameters:**
- `confirm`: Whether to ask for user confirmation

## Error Handling

The class includes comprehensive error handling for:
- File not found errors
- Missing CSV columns
- Database connection issues  
- Data validation errors
- Duplicate handling

All methods return structured results with success indicators and descriptive error messages.

## Integration with TidyScreen

ChemSpace seamlessly integrates with the TidyScreen project management system:

1. **Project Object Integration**: Uses ActivateProject objects directly for project access
2. **Attribute Access**: Accesses all project attributes (name, path, description, etc.) from the project object
3. **Path Management**: Uses project path from the ActivateProject object for database location
4. **Project Validation**: Ensures project exists through the ActivateProject object validation
5. **Directory Structure**: Follows TidyScreen project conventions automatically

## Best Practices

1. **Project Setup**: Always create the project using `CreateProject` before activating it with `ActivateProject`
2. **Object Management**: Use `ActivateProject` objects to initialize `ChemSpace` for better integration
3. **CSV Format**: Ensure CSV files have proper headers matching the expected column names
4. **Data Validation**: Review import results and handle errors appropriately
5. **Backup**: Consider exporting data periodically for backup purposes
6. **Performance**: Use filters and limits for large datasets
7. **Error Handling**: Always check if the project exists before creating ChemSpace instances

## Example CSV Format

```csv
name,smiles,flag
aspirin,CC(=O)OC1=CC=CC=C1C(=O)O,approved
caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C,approved
experimental_1,CCCN(CC)CC(=O)NC1=CC=CC=C1,experimental
```

## Dependencies

- `pandas`: For CSV file processing
- `sqlite3`: For database operations (built-in)
- `tidyscreen`: For project management integration
- Standard library modules: `os`, `csv`, `datetime`, `typing`
