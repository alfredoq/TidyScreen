# TidyScreen Package - UML Class Diagram Documentation

## Package Overview
The TidyScreen package is organized into several modules that work together to manage project creation, activation, and deletion with database persistence and file system operations.

## Class Hierarchy

### 1. DatabaseManager (Base Class)
**Location:** `src/databases/DatabaseManager.py`
**Purpose:** Provides base functionality for SQLite database operations

**Attributes:**
- None (stateless base class)

**Methods:**
- `+ check_db(db: str) : bool` - Verifies database existence
- `+ create_projects_database(db_path: str) : bool` - Creates project database
- `+ connect_db(db: str) : sqlite3.Connection` - Establishes database connection
- `+ execute_query(conn, query: str, message: str) : void` - Executes SQL queries

### 2. ProjectsManagement (Inherits from DatabaseManager)
**Location:** `src/projects/projects_management.py`
**Purpose:** Manages all project-related database operations (CRUD)

**Attributes:**
- `- __projects_db: str` - Path to projects database

**Methods:**
- `+ __init__()` - Initializes with database path
- `+ list_all_projects(print_output: bool = True) : list | None` - Lists all projects
- `+ load_project_info(project_name: str) : dict` - Retrieves project information
- `+ create_project_entry(name: str, path: str, description: str = None) : dict` - Creates database entry
- `+ delete_project_entry(project_name: str, delete_directory: bool = True) : dict` - Deletes project

### 3. ActivateProject (Inherits from DatabaseManager)
**Location:** `src/tidyscreen.py`
**Purpose:** Loads and manages existing projects from database

**Attributes:**
- `- __projects_db: str` - Path to projects database
- `+ name: str` - Project name
- `+ id: int` - Project database ID
- `+ path: str` - Project file system path
- `+ description: str` - Project description
- `+ created_date: str` - Project creation timestamp
- `+ _project_exists: bool` - Project existence flag

**Methods:**
- `+ __init__(name: str)` - Initializes and loads project
- `+ load_project_info() : void` - Loads project data from database
- `+ get_project_attributes() : dict` - Returns all project attributes
- `+ project_exists() : bool` - Checks if project exists

### 4. CreateProject (Inherits from DatabaseManager)
**Location:** `src/tidyscreen.py`
**Purpose:** Creates new projects with directory structure and database entry

**Attributes:**
- `- __projects_db: str` - Path to projects database
- `+ name: str` - Project name
- `+ path: str` - Project file system path
- `+ description: str` - Project description
- `+ id: int` - Project database ID (after creation)
- `+ created_date: str` - Project creation timestamp
- `+ _directory_created: bool` - Directory creation success flag
- `+ _project_created: bool` - Database entry creation success flag

**Methods:**
- `+ __init__(name: str, path: str, description: str = None)` - Creates new project
- `+ load_project_template() : dict` - Loads JSON template configuration
- `+ _get_default_structure() : dict` - Provides fallback structure
- `+ create_project_directory() : void` - Creates directory structure
- `+ _create_folders(folders: list) : void` - Creates main folders
- `+ _create_subfolders(subfolders: dict) : void` - Creates nested folders
- `+ _create_template_files(files: dict) : void` - Creates template files
- `+ create_project_entry() : void` - Creates database entry
- `+ directory_created_successfully() : bool` - Checks creation status

### 5. DeleteProject (Inherits from DatabaseManager)
**Location:** `src/tidyscreen.py`
**Purpose:** Deletes projects from database and optionally from file system

**Attributes:**
- `- __projects_db: str` - Path to projects database
- `+ name: str` - Project name to delete
- `+ delete_directory: bool` - Whether to delete file system directory
- `+ confirm: bool` - Whether to show confirmation prompt
- `+ id: int` - Project database ID
- `+ path: str` - Project file system path
- `+ description: str` - Project description
- `+ created_date: str` - Project creation timestamp
- `+ directory_deleted: bool` - Directory deletion success flag
- `+ _project_deleted: bool` - Project deletion success flag

**Methods:**
- `+ __init__(name: str, delete_directory: bool = True, confirm: bool = True)` - Initializes deletion
- `+ delete_project_entry() : void` - Performs project deletion with confirmation
- `+ project_deleted_successfully() : bool` - Checks deletion status

### 6. ProjectTemplate (Configuration)
**Location:** `src/config/project_template.json`
**Purpose:** JSON configuration defining project directory structure

**Structure:**
```json
{
  "project_structure": {
    "folders": ["src", "tests", "docs", ...],
    "subfolders": {
      "src": ["models", "views", "controllers", ...],
      ...
    },
    "files": {
      "README.md": "template content...",
      ...
    }
  }
}
```

### 7. projects() Function
**Location:** `src/tidyscreen.py`
**Purpose:** Standalone function to list all projects

**Signature:** `+ projects() : void`
**Behavior:** Creates ProjectsManagement instance and calls list_all_projects()

## Relationships

### Inheritance
- `ProjectsManagement` → `DatabaseManager`
- `ActivateProject` → `DatabaseManager`
- `CreateProject` → `DatabaseManager`
- `DeleteProject` → `DatabaseManager`

### Composition/Usage
- `ActivateProject` uses `ProjectsManagement` for data loading
- `CreateProject` uses `ProjectsManagement` for database operations
- `DeleteProject` uses `ProjectsManagement` for database operations
- `CreateProject` reads `ProjectTemplate` JSON configuration
- `projects()` function uses `ProjectsManagement`

### Dependencies
- All classes depend on SQLite database for persistence
- `CreateProject` and `DeleteProject` interact with file system
- `CreateProject` depends on JSON configuration file

## Design Patterns

1. **Inheritance:** All main classes inherit from DatabaseManager for common database functionality
2. **Composition:** Classes use ProjectsManagement for specialized database operations
3. **Template Method:** CreateProject uses JSON templates for flexible project structure
4. **Strategy:** Different project operations (create, activate, delete) are encapsulated in separate classes
5. **Factory-like:** projects() function provides a simple interface to project listing

## Key Features

1. **Database Abstraction:** DatabaseManager provides consistent SQLite interface
2. **Flexible Project Structure:** JSON-based templates allow customizable directory hierarchies
3. **Error Handling:** Comprehensive error handling and user feedback
4. **Safety Features:** Confirmation prompts for destructive operations
5. **Separation of Concerns:** Clear separation between database operations, file system operations, and business logic
