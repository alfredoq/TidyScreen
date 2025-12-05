classDiagram
    %% TidyScreen Package Class Diagram

    %% Base Class
    class DatabaseManager {
        +check_db(db: str) bool
        +create_projects_database(db_path: str) bool
        +connect_db(db: str) Connection
        +execute_query(conn, query: str, message: str) void
    }

    %% Projects Management
    class ProjectsManagement {
        -__projects_db: str
        +__init__()
        +list_all_projects(print_output: bool) list|None
        +load_project_info(project_name: str) dict
        +create_project_entry(name: str, path: str, description: str) dict
        +delete_project_entry(project_name: str, delete_directory: bool) dict
    }

    %% Main Classes
    class ActivateProject {
        -__projects_db: str
        +name: str
        +id: int
        +path: str
        +description: str
        +created_date: str
        +_project_exists: bool
        +__init__(name: str)
        +load_project_info() void
        +get_project_attributes() dict
        +project_exists() bool
    }

    class CreateProject {
        -__projects_db: str
        +name: str
        +path: str
        +description: str
        +id: int
        +created_date: str
        +_directory_created: bool
        +_project_created: bool
        +__init__(name: str, path: str, description: str)
        +load_project_template() dict
        +_get_default_structure() dict
        +create_project_directory() void
        +_create_folders(folders: list) void
        +_create_subfolders(subfolders: dict) void
        +_create_template_files(files: dict) void
        +create_project_entry() void
        +directory_created_successfully() bool
    }

    class DeleteProject {
        -__projects_db: str
        +name: str
        +delete_directory: bool
        +confirm: bool
        +id: int
        +path: str
        +description: str
        +created_date: str
        +directory_deleted: bool
        +_project_deleted: bool
        +__init__(name: str, delete_directory: bool, confirm: bool)
        +delete_project_entry() void
        +project_deleted_successfully() bool
    }

    %% Configuration
    class ProjectTemplate {
        <<JSON>>
        +project_structure: object
        +folders: array
        +subfolders: object
        +files: object
    }

    %% Function
    class projects_function["projects()"] {
        <<function>>
        +projects() void
    }

    %% Inheritance
    DatabaseManager <|-- ProjectsManagement
    DatabaseManager <|-- ActivateProject
    DatabaseManager <|-- CreateProject
    DatabaseManager <|-- DeleteProject

    %% Composition/Usage
    ActivateProject ..> ProjectsManagement : uses
    CreateProject ..> ProjectsManagement : uses
    DeleteProject ..> ProjectsManagement : uses
    CreateProject ..> ProjectTemplate : reads
    projects_function ..> ProjectsManagement : uses

    %% Notes
    note for DatabaseManager "Base class for SQLite database operations"
    note for ProjectsManagement "Handles project CRUD operations"
    note for CreateProject "Creates projects with JSON-based directory structure"
    note for DeleteProject "Deletes projects with confirmation"
    note for ActivateProject "Loads and manages existing projects"
    note for ProjectTemplate "JSON configuration for project templates"
