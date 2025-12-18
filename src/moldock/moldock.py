import os
from tidyscreen import tidyscreen
import sys
from typing import Dict, List, Tuple, Optional, Any
from rdkit.Chem import AllChem
from rdkit import Chem

# Try to import tqdm for progress bars
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

# Add the parent directory to path to import our local tidyscreen module
parent_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, parent_dir)

# Import from our local tidyscreen module
from tidyscreen import tidyscreen
ActivateProject = tidyscreen.ActivateProject

class MolDock:
    """
    MolDock class for managing molecular docking assays within a project.
    Uses an ActivateProject object to access project information and database functionality.
    
    """
    
    def __init__(self, project_obj: ActivateProject):
        """
        Initialize MolDock with an ActivateProject object.
        
        Args:
            project_obj (ActivateProject): An instantiated ActivateProject object
        """
        # Validate that we received a proper ActivateProject object
        if not isinstance(project_obj, ActivateProject):
            raise TypeError("ChemSpace requires an ActivateProject object")
        
        # Check if project exists
        if not project_obj.project_exists():
            raise ValueError(f"Project '{project_obj.name}' not found. Please create the project first.")
        
        # Store the project object and its attributes
        self.project = project_obj
        self.name = project_obj.name
        self.path = project_obj.path
        self.description = getattr(project_obj, 'description', None)
        self.id = getattr(project_obj, 'id', None)
        self.created_date = getattr(project_obj, 'created_date', None)
        
        # Set up chemspace database path within the project directory
        self.__chemspace_db = os.path.join(self.path, 'chemspace/processed_data', 'chemspace.db')
        
        # Set up docking assays registers database path within the project directory
        self.__docking_registers_db = os.path.join(self.path, 'docking/docking_registers', 'registers.db')

        # Set up docking params registers database path within the project directory
        self.__docking_params_db = os.path.join(self.path, 'docking/params', 'params.db')

        # Ensure data directory exists
        data_dir = os.path.dirname(self.__chemspace_db)
        os.makedirs(data_dir, exist_ok=True)

    def dock_table(self, show_details: bool = True, 
                filter_by_type: Optional[str] = None,
                sort_by: str = 'name',
                show_all_tables: bool = False) -> Optional[Any]:
        """
        List tables available in the chemspace database and allow selection for docking preparation.
        By default, only shows tables ready for docking (with SDF blobs).
        
        Args:
            show_details (bool): Whether to show detailed information including compound counts
            filter_by_type (Optional[str]): Filter tables by type ('original', 'reaction_products', 'filtered_compounds', etc.)
            sort_by (str): Sort criteria ('name', 'compounds', 'type', 'date')
            show_all_tables (bool): If True, shows all tables regardless of docking readiness
            
        Returns:
            Optional[Any]: Selected table name, list of SMILES, or None if cancelled
        """
        try:
            print(f"🧬 MOLDOCK - Table Selection for Docking Preparation")
            print("=" * 70)
            
            # Check if chemspace database exists
            if not os.path.exists(self.__chemspace_db):
                print(f"❌ ChemSpace database not found: {self.__chemspace_db}")
                print("   Please ensure ChemSpace has been initialized for this project.")
                return None
            
            # Connect to chemspace database and get available tables
            try:
                import sqlite3
                conn = sqlite3.connect(self.__chemspace_db)
                cursor = conn.cursor()
                
                # Get all tables (reusing ChemSpace pattern)
                cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")
                available_tables = [row[0] for row in cursor.fetchall()]
                conn.close()
                
            except Exception as e:
                print(f"❌ Error accessing ChemSpace database: {e}")
                return None
            
            if not available_tables:
                print("📋 No tables found in ChemSpace database")
                print("   Load some compound data into ChemSpace first.")
                return None
            
            print(f"📊 Found {len(available_tables)} table(s) in ChemSpace database")
            
            # Collect table information using existing patterns
            table_data = []
            ready_table_data = []
            not_ready_count = 0
            
            if show_details:
                print("📊 Analyzing tables for docking preparation...")
            
            for table_name in available_tables:
                try:
                    # Get compound count (reusing ChemSpace pattern)
                    table_info = self._get_table_info_for_docking(table_name)
                    
                    # Apply type filter if specified
                    if filter_by_type and table_info['type'] != filter_by_type:
                        continue
                    
                    # Separate ready vs not ready tables
                    if table_info['docking_status'] == 'ready for docking':
                        ready_table_data.append(table_info)
                    else:
                        not_ready_count += 1
                    
                    # Add to all tables list (for show_all_tables option)
                    table_data.append(table_info)
                    
                except Exception as e:
                    print(f"   ⚠️  Error analyzing table '{table_name}': {e}")
                    continue
            
            # Show summary of table readiness
            total_analyzed = len(table_data)
            ready_count = len(ready_table_data)
            
            print(f"\n📋 TABLE READINESS SUMMARY:")
            print(f"   📊 Total tables analyzed: {total_analyzed}")
            print(f"   ✅ Ready for docking: {ready_count}")
            print(f"   ❌ Not ready for docking: {not_ready_count}")
            
            # Determine which tables to display
            if not show_all_tables:
                # Default behavior: show only ready tables
                tables_to_display = ready_table_data
                
                if not ready_table_data:
                    print(f"\n⚠️  NO TABLES READY FOR DOCKING")
                    print(f"   All {total_analyzed} tables need 3D molecular preparation")
                    print(f"   💡 Suggestion: Use ChemSpace.generate_mols_in_table() to prepare molecules")
                    
                    # Offer to show all tables
                    if total_analyzed > 0:
                        print(f"\n🔄 OPTIONS:")
                        print(f"   1. Show all tables (including not ready)")
                        print(f"   2. Cancel and prepare molecules first")
                        
                        while True:
                            try:
                                choice = input(f"\n🧬 Choose option (1/2): ").strip()
                                
                                if choice == '1':
                                    print(f"\n📋 Showing all {total_analyzed} tables (including not ready)...")
                                    tables_to_display = table_data
                                    show_all_tables = True
                                    break
                                elif choice == '2' or choice.lower() in ['cancel', 'quit', 'exit']:
                                    print("❌ Docking preparation cancelled")
                                    return None
                                else:
                                    print("❌ Invalid choice. Please enter 1 or 2")
                                    continue
                                    
                            except KeyboardInterrupt:
                                print("\n❌ Docking preparation cancelled")
                                return None
                    else:
                        return None
                else:
                    print(f"\n✅ Showing {ready_count} table(s) ready for docking")
                    
            else:
                # Show all tables when explicitly requested
                tables_to_display = table_data
                print(f"\n📋 Showing all {total_analyzed} tables (ready + not ready)")
            
            if not tables_to_display:
                if filter_by_type:
                    print(f"❌ No tables found with type '{filter_by_type}'")
                else:
                    print("❌ No valid tables found for docking preparation")
                return None
            
            # Sort tables (reusing existing logic)
            tables_to_display = self._sort_table_data_for_docking(tables_to_display, sort_by)
            
            # Display table list for selection with readiness context
            selected_table = self._display_and_select_docking_table(
                tables_to_display, filter_by_type, show_all_tables, ready_count, not_ready_count
            )
            
            if selected_table:
                print(f"\n✅ Selected table for docking: '{selected_table}'")
                
                # Show preparation summary and get user choice
                user_choice = self._show_docking_preparation_summary(selected_table, tables_to_display)
                
                if user_choice == 1:
                    # Return table name only
                    return selected_table
                elif user_choice == 2:
                    # Return SMILES list for immediate processing
                    print(f"\n🧪 Loading SMILES from table '{selected_table}'...")
                    smiles_data = self._get_smiles_from_table(selected_table)
                    
                    if smiles_data:
                        print(f"✅ Loaded {len(smiles_data)} compounds with SMILES")
                        return smiles_data
                    else:
                        print("❌ No valid SMILES found in selected table")
                        return None
                else:
                    return None
            else:
                print("❌ No table selected for docking preparation")
                return None
            
        except Exception as e:
            print(f"❌ Error in dock_table method: {e}")
            return None

    def _show_docking_preparation_summary(self, selected_table: str, table_data: List[Dict]) -> Optional[int]:
        """
        Show summary of selected table for docking preparation and get user choice.
        
        Args:
            selected_table (str): Name of selected table
            table_data (List[Dict]): Table data list
            
        Returns:
            Optional[int]: User choice (1=table name, 2=SMILES list, None=cancel)
        """
        try:
            # Find selected table info
            table_info = next((t for t in table_data if t['name'] == selected_table), None)
            
            if not table_info:
                return None
            
            print(f"\n🎯 DOCKING PREPARATION SUMMARY")
            print("-" * 50)
            print(f"📋 Selected Table: '{selected_table}'")
            print(f"📊 Compounds: {table_info['compound_count']:,}")
            print(f"🎯 Readiness: {table_info['docking_readiness']['status']}")
            
            readiness = table_info['docking_readiness']
            if readiness['score'] >= 70:
                print("✅ Table is ready for docking preparation")
            elif readiness['score'] >= 50:
                print("🟡 Table has minor issues but can be used")
            else:
                print("⚠️  Table has significant issues for docking")
            
            print(f"\n🔄 PREPARATION OPTIONS:")
            print("   1. Return table name (for later processing)")
            print("   2. Load SMILES immediately (for molecule preparation)")
            print("   3. Cancel selection")
            print("-" * 50)
            
            # Get user choice
            while True:
                try:
                    choice = input("🧬 Select preparation option (1/2/3): ").strip()
                    
                    if choice == '1':
                        print("✅ Will return table name for later processing")
                        return 1
                    elif choice == '2':
                        print("🧪 Will load SMILES for immediate molecule preparation")
                        return 2
                    elif choice == '3' or choice.lower() in ['cancel', 'quit', 'exit']:
                        print("❌ Preparation cancelled")
                        return None
                    else:
                        print("❌ Invalid choice. Please enter 1, 2, or 3")
                        continue
                        
                except KeyboardInterrupt:
                    print("\n❌ Preparation cancelled")
                    return None
            
        except Exception as e:
            print(f"⚠️  Error showing preparation summary: {e}")
            return None

    def _get_smiles_from_table(self, table_name: str) -> Optional[List[Dict[str, Any]]]:
        """
        Extract SMILES and associated data from the selected table.
        
        Args:
            table_name (str): Name of the table to extract SMILES from
            
        Returns:
            Optional[List[Dict]]: List of compound dictionaries with SMILES and metadata
        """
        try:
            import sqlite3
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Get table schema to identify available columns
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns = [row[1] for row in cursor.fetchall()]
            
            # Check for required and optional columns
            has_smiles = 'smiles' in [col.lower() for col in columns]
            if not has_smiles:
                print("❌ No SMILES column found in table")
                conn.close()
                return None
            
            # Build query to get available data
            base_columns = ['smiles']
            optional_columns = ['name', 'id', 'inchi_key', 'molecular_weight', 'logp']
            
            available_columns = []
            for col in columns:
                col_lower = col.lower()
                if col_lower == 'smiles':
                    available_columns.append('smiles')
                elif col_lower in ['name', 'compound_name', 'title']:
                    available_columns.append(f"{col} as name")
                elif col_lower in ['id', 'compound_id', 'mol_id']:
                    available_columns.append(f"{col} as id")
                elif col_lower in ['inchi_key', 'inchikey', 'inchi']:
                    available_columns.append(f"{col} as inchi_key")
                elif col_lower in ['molecular_weight', 'mw', 'mol_weight']:
                    available_columns.append(f"{col} as molecular_weight")
                elif col_lower in ['logp', 'clogp', 'alogp']:
                    available_columns.append(f"{col} as logp")
            
            # Execute query
            query = f"SELECT {', '.join(available_columns)} FROM {table_name} WHERE smiles IS NOT NULL AND smiles != ''"
            
            if TQDM_AVAILABLE:
                # First get count for progress bar
                cursor.execute(f"SELECT COUNT(*) FROM {table_name} WHERE smiles IS NOT NULL AND smiles != ''")
                total_count = cursor.fetchone()[0]
                
                cursor.execute(query)
                results = cursor.fetchall()
                
                progress_bar = tqdm(total=len(results), desc="Loading SMILES", unit="compounds")
            else:
                cursor.execute(query)
                results = cursor.fetchall()
                progress_bar = None
                print(f"📊 Processing {len(results)} compounds...")
            
            smiles_data = []
            processed_count = 0
            
            for row in results:
                try:
                    # Create compound dictionary
                    compound_dict = {}
                    
                    for i, col_info in enumerate(available_columns):
                        if ' as ' in col_info:
                            col_name = col_info.split(' as ')[1]
                        else:
                            col_name = col_info
                        
                        compound_dict[col_name] = row[i] if i < len(row) else None
                    
                    # Validate SMILES
                    smiles = compound_dict.get('smiles', '').strip()
                    if smiles and len(smiles) > 0:
                        # Add metadata for docking preparation
                        compound_dict['conf_rank'] = 0  # Default conformer rank
                        compound_dict['source_table'] = table_name
                        compound_dict['compound_index'] = processed_count
                        
                        smiles_data.append(compound_dict)
                        processed_count += 1
                    
                    if progress_bar:
                        progress_bar.update(1)
                    elif processed_count % 1000 == 0:
                        print(f"   📊 Processed: {processed_count:,} compounds")
                    
                except Exception as e:
                    print(f"   ⚠️  Error processing compound at index {len(smiles_data)}: {e}")
                    continue
            
            if progress_bar:
                progress_bar.close()
            
            conn.close()
            
            print(f"\n✅ Successfully loaded {len(smiles_data)} valid compounds")
            print(f"   📋 Table: {table_name}")
            print(f"   🧪 Valid SMILES: {len(smiles_data):,}")
            
            if len(smiles_data) == 0:
                print("❌ No valid SMILES found in table")
                return None
            
            return smiles_data
            
        except Exception as e:
            print(f"❌ Error loading SMILES from table: {e}")
            return None

    def _get_table_info_for_docking(self, table_name: str) -> Dict[str, Any]:
        """
        Get table information specifically for docking preparation.
        Reuses existing database access patterns.
        
        Args:
            table_name (str): Name of the table to analyze
            
        Returns:
            Dict[str, Any]: Table information dictionary
        """
        try:
            import sqlite3
            conn = sqlite3.connect(self.__chemspace_db)
            cursor = conn.cursor()
            
            # Get basic table info
            cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
            compound_count = cursor.fetchone()[0]
            
            # Check table schema for required columns
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns = [row[1] for row in cursor.fetchall()]
            has_smiles = 'smiles' in [col.lower() for col in columns]
            has_sdf_blob = 'sdf_blob' in [col.lower() for col in columns]
            
            # Count SDF blobs if column exists
            sdf_blob_count = 0
            if has_sdf_blob:
                cursor.execute(f"SELECT COUNT(*) FROM {table_name} WHERE sdf_blob IS NOT NULL")
                sdf_blob_count = cursor.fetchone()[0]
            
            # Get sample data to assess docking suitability
            if has_smiles and compound_count > 0:
                cursor.execute(f"SELECT smiles FROM {table_name} WHERE smiles IS NOT NULL LIMIT 5")
                sample_smiles = [row[0] for row in cursor.fetchall() if row[0] and row[0].strip()]
                valid_smiles_count = len(sample_smiles)
            else:
                valid_smiles_count = 0
                sample_smiles = []
            
            conn.close()
            
            # Classify table type (reusing existing logic from ChemSpace)
            table_type = self._classify_table_type_for_docking(table_name)
            
            # Assess docking readiness with SDF blob information
            docking_readiness = self._assess_docking_readiness(
                compound_count, has_smiles, valid_smiles_count, sample_smiles,
                has_sdf_blob, sdf_blob_count
            )
            
            return {
                'name': table_name,
                'type': table_type,
                'compound_count': compound_count,
                'has_smiles': has_smiles,
                'has_sdf_blob': has_sdf_blob,
                'sdf_blob_count': sdf_blob_count,
                'valid_smiles_count': valid_smiles_count,
                'docking_readiness': docking_readiness,
                'readiness_score': docking_readiness['score'],
                'docking_status': docking_readiness['docking_status'],
                'icon': self._get_table_icon_for_docking(table_type, docking_readiness['score'], has_sdf_blob, sdf_blob_count)
            }
            
        except Exception as e:
            return {
                'name': table_name,
                'type': 'error',
                'compound_count': 0,
                'has_smiles': False,
                'has_sdf_blob': False,
                'sdf_blob_count': 0,
                'valid_smiles_count': 0,
                'docking_readiness': {'score': 0, 'status': 'error', 'docking_status': 'NOT ready for docking', 'issues': [str(e)]},
                'readiness_score': 0,
                'docking_status': 'NOT ready for docking',
                'icon': '❌'
            }

    def _classify_table_type_for_docking(self, table_name: str) -> str:
        """
        Classify table type for docking purposes.
        Reuses existing classification logic from ChemSpace.
        
        Args:
            table_name (str): Name of the table
            
        Returns:
            str: Table type classification
        """
        name_lower = table_name.lower()
        
        # Reuse ChemSpace classification patterns
        if '_products_' in name_lower or 'reaction' in name_lower:
            return 'reaction_products'
        elif '_filtered_' in name_lower or 'filter' in name_lower:
            return 'filtered_compounds'
        elif 'stereoisomer' in name_lower or 'stereo' in name_lower:
            return 'stereoisomers'
        elif 'ersilia' in name_lower or 'prediction' in name_lower:
            return 'predicted_compounds'
        elif name_lower.startswith('temp_') or name_lower.startswith('tmp_'):
            return 'temporary'
        elif any(word in name_lower for word in ['workflow', 'batch', 'processed']):
            return 'processed'
        else:
            return 'original'

    def _assess_docking_readiness(self, compound_count: int, has_smiles: bool, 
                                valid_smiles_count: int, sample_smiles: List[str],
                                has_sdf_blob: bool = False, sdf_blob_count: int = 0) -> Dict[str, Any]:
        """
        Assess how ready a table is for docking preparation.
        
        Args:
            compound_count (int): Total number of compounds
            has_smiles (bool): Whether table has SMILES column
            valid_smiles_count (int): Number of valid SMILES in sample
            sample_smiles (List[str]): Sample SMILES strings
            has_sdf_blob (bool): Whether table has SDF blob column
            sdf_blob_count (int): Number of compounds with SDF blobs
            
        Returns:
            Dict[str, Any]: Readiness assessment
        """
        try:
            issues = []
            score = 100  # Start with perfect score and deduct points
            docking_status = "NOT ready for docking"  # Default status
            
            # Check basic requirements
            if compound_count == 0:
                issues.append("Table is empty")
                score = 0
            elif compound_count < 5:
                issues.append(f"Very few compounds ({compound_count})")
                score -= 20
            
            if not has_smiles:
                issues.append("No SMILES column found")
                score = 0
            elif valid_smiles_count == 0:
                issues.append("No valid SMILES found in sample")
                score -= 50
            elif valid_smiles_count < len(sample_smiles):
                issues.append("Some invalid SMILES detected")
                score -= 10
            
            # Check for SDF blob column (critical for docking readiness)
            if has_sdf_blob:
                if sdf_blob_count > 0:
                    docking_status = "ready for docking"
                    print(f"   ✅ Found {sdf_blob_count:,} prepared molecules with 3D coordinates")
                else:
                    issues.append("SDF blob column exists but no prepared molecules found")
                    docking_status = "NOT ready for docking"
                    score -= 30
            else:
                issues.append("No SDF blob column - molecules need 3D preparation")
                docking_status = "NOT ready for docking"
                score -= 40
            
            # Check SMILES complexity for docking suitability (only if we have valid SMILES)
            if sample_smiles and score > 0:
                complexity_assessment = self._assess_smiles_complexity(sample_smiles)
                if complexity_assessment['has_issues']:
                    issues.extend(complexity_assessment['issues'])
                    score -= complexity_assessment['score_penalty']
            
            # Determine overall status with docking-specific criteria
            if has_sdf_blob and sdf_blob_count > 0:
                if score >= 90:
                    status = 'excellent'
                elif score >= 70:
                    status = 'good'
                elif score >= 50:
                    status = 'fair'
                else:
                    status = 'poor'
            else:
                # Without SDF blobs, maximum status is 'needs_preparation'
                if score >= 70:
                    status = 'needs_preparation'
                elif score >= 50:
                    status = 'fair'
                elif score > 0:
                    status = 'poor'
                else:
                    status = 'unsuitable'
            
            return {
                'score': max(0, score),
                'status': status,
                'docking_status': docking_status,
                'issues': issues,
                'compound_count': compound_count,
                'has_smiles': has_smiles,
                'has_sdf_blob': has_sdf_blob,
                'sdf_blob_count': sdf_blob_count,
                'valid_smiles_sample': valid_smiles_count
            }
            
        except Exception as e:
            return {
                'score': 0,
                'status': 'error',
                'docking_status': 'NOT ready for docking',
                'issues': [f"Assessment error: {e}"],
                'compound_count': compound_count,
                'has_smiles': has_smiles,
                'has_sdf_blob': False,
                'sdf_blob_count': 0,
                'valid_smiles_sample': 0
            }

    def _assess_smiles_complexity(self, sample_smiles: List[str]) -> Dict[str, Any]:
        """
        Assess SMILES complexity for docking suitability.
        
        Args:
            sample_smiles (List[str]): Sample SMILES strings
            
        Returns:
            Dict[str, Any]: Complexity assessment
        """
        try:
            issues = []
            score_penalty = 0
            
            for smiles in sample_smiles[:3]:  # Check first 3 samples
                # Check for very simple molecules (might not be drug-like)
                if len(smiles) < 10:
                    issues.append("Very simple molecules detected")
                    score_penalty += 5
                
                # Check for very complex molecules (might be problematic for docking)
                elif len(smiles) > 200:
                    issues.append("Very complex molecules detected")
                    score_penalty += 10
                
                # Check for inorganic/problematic patterns
                problematic_patterns = ['[Na+]', '[K+]', '[Cl-]', '[SO4-2]', '[PO4-3]']
                if any(pattern in smiles for pattern in problematic_patterns):
                    issues.append("Inorganic/salt components detected")
                    score_penalty += 15
            
            return {
                'has_issues': len(issues) > 0,
                'issues': list(set(issues)),  # Remove duplicates
                'score_penalty': min(score_penalty, 30)  # Cap penalty
            }
            
        except Exception:
            return {
                'has_issues': True,
                'issues': ["Error assessing SMILES complexity"],
                'score_penalty': 10
            }

    def _get_table_icon_for_docking(self, table_type: str, readiness_score: int, 
                                has_sdf_blob: bool = False, sdf_blob_count: int = 0) -> str:
        """
        Get appropriate icon for table based on type, docking readiness, and SDF blob availability.
        
        Args:
            table_type (str): Table type
            readiness_score (int): Docking readiness score
            has_sdf_blob (bool): Whether table has SDF blob column
            sdf_blob_count (int): Number of compounds with SDF blobs
            
        Returns:
            str: Emoji icon
        """
        # Base icons by type (reusing ChemSpace patterns)
        type_icons = {
            'original': '📋',
            'reaction_products': '🧪',
            'filtered_compounds': '🔍',
            'stereoisomers': '🔄',
            'predicted_compounds': '🧠',
            'processed': '⚗️',
            'temporary': '📄',
            'error': '❌'
        }
        
        base_icon = type_icons.get(table_type, '📊')
        
        # Docking readiness indicators with SDF blob consideration
        if has_sdf_blob and sdf_blob_count > 0:
            # Has prepared 3D molecules - ready for docking
            if readiness_score >= 90:
                return f"{base_icon}🚀"  # Ready to dock - excellent
            elif readiness_score >= 70:
                return f"{base_icon}✅"  # Ready to dock - good
            else:
                return f"{base_icon}🟢"  # Ready to dock - fair quality
        else:
            # No 3D molecules prepared - needs preparation
            if readiness_score >= 70:
                return f"{base_icon}🔧"  # Needs 3D preparation
            elif readiness_score >= 50:
                return f"{base_icon}🟡"  # Needs preparation + has issues
            elif readiness_score > 0:
                return f"{base_icon}🟠"  # Poor quality, needs work
            else:
                return f"{base_icon}🔴"  # Unsuitable for docking

    def _sort_table_data_for_docking(self, table_data: List[Dict], sort_by: str) -> List[Dict]:
        """
        Sort table data for docking display.
        Reuses existing sorting logic with docking-specific priorities.
        
        Args:
            table_data (List[Dict]): List of table information dictionaries
            sort_by (str): Sort criteria
            
        Returns:
            List[Dict]: Sorted table data
        """
        try:
            if sort_by == 'readiness':
                # Sort by docking readiness score (highest first)
                return sorted(table_data, key=lambda x: x['readiness_score'], reverse=True)
            elif sort_by == 'compounds':
                return sorted(table_data, key=lambda x: x['compound_count'], reverse=True)
            elif sort_by == 'type':
                return sorted(table_data, key=lambda x: (x['type'], x['name'].lower()))
            else:  # sort by name (default)
                return sorted(table_data, key=lambda x: x['name'].lower())
        except Exception:
            return table_data

    def _display_and_select_docking_table(self, table_data: List[Dict], 
                                        filter_type: Optional[str],
                                        show_all_tables: bool = False,
                                        ready_count: int = 0,
                                        not_ready_count: int = 0) -> Optional[str]:
        """
        Display table list and handle user selection for docking.
        
        Args:
            table_data (List[Dict]): Table information data
            filter_type (Optional[str]): Applied filter type
            show_all_tables (bool): Whether showing all tables or only ready ones
            ready_count (int): Total number of ready tables
            not_ready_count (int): Total number of not ready tables
            
        Returns:
            Optional[str]: Selected table name or None if cancelled
        """
        try:
            # Enhanced header showing filtering context
            if show_all_tables:
                print(f"📊 Displaying ALL tables: {len(table_data)} (Ready: {ready_count}, Not Ready: {not_ready_count})")
                if filter_type:
                    print(f"🔍 Additionally filtered by type: '{filter_type}'")
            else:
                print(f"✅ Displaying READY tables only: {len(table_data)} of {ready_count + not_ready_count} total")
                if filter_type:
                    print(f"🔍 Additionally filtered by type: '{filter_type}'")
            
            print("-" * 90)
            print(f"{'#':<3} {'Icon':<6} {'Table Name':<25} {'Type':<15} {'Compounds':<10} {'3D Mols':<8} {'Docking Status':<15}")
            print("-" * 90)
            
            for i, table_info in enumerate(table_data, 1):
                compounds_str = f"{table_info['compound_count']:,}"
                sdf_count_str = f"{table_info['sdf_blob_count']:,}" if table_info['has_sdf_blob'] else "0"
                docking_status = "✅ Ready" if table_info['docking_status'] == 'ready for docking' else "❌ Not Ready"
                
                print(f"{i:<3} {table_info['icon']:<6} {table_info['name'][:24]:<25} "
                    f"{table_info['type'][:14]:<15} {compounds_str:<10} {sdf_count_str:<8} {docking_status:<15}")
            
            print("-" * 90)
            
            # Show enhanced readiness legend
            print("📊 Icons: 🚀✅Ready to dock  🔧Needs 3D prep  🟡🟠Issues  🔴Unsuitable")
            print("💡 3D Mols: Number of molecules with prepared 3D coordinates (SDF blobs)")
            
            # Show filtering status and options
            if not show_all_tables and not_ready_count > 0:
                print(f"ℹ️  Note: {not_ready_count} not-ready tables are hidden. Use show_all_tables=True to see them.")
            
            print("\nCommands: Enter table number, table name, 'details <num>' for info, or 'cancel'")
            
            # Rest of the selection logic remains the same...
            while True:
                try:
                    selection = input(f"\n🧬 Select table for docking preparation: ").strip()
                    
                    if selection.lower() in ['cancel', 'quit', 'exit']:
                        return None
                    
                    # Handle details command
                    if selection.lower().startswith('details '):
                        try:
                            detail_idx = int(selection.split()[1]) - 1
                            if 0 <= detail_idx < len(table_data):
                                self._show_table_details_for_docking(table_data[detail_idx])
                                continue
                            else:
                                print(f"❌ Invalid table number. Please enter 1-{len(table_data)}")
                                continue
                        except (IndexError, ValueError):
                            print("❌ Invalid details command. Use 'details <number>'")
                            continue
                    
                    # Try as number first
                    try:
                        table_idx = int(selection) - 1
                        if 0 <= table_idx < len(table_data):
                            selected_info = table_data[table_idx]
                            
                            # Enhanced warning for tables not ready for docking
                            if selected_info['docking_status'] == 'NOT ready for docking':
                                print(f"\n⚠️  WARNING: Table '{selected_info['name']}' is NOT ready for docking")
                                if not selected_info['has_sdf_blob']:
                                    print("   🔧 Missing: 3D molecular coordinates (SDF blobs)")
                                    print("   💡 Suggestion: Use ChemSpace.generate_mols_in_table() to prepare 3D molecules")
                                elif selected_info['sdf_blob_count'] == 0:
                                    print("   📊 SDF blob column exists but no prepared molecules found")
                                
                                print("Issues:")
                                for issue in selected_info['docking_readiness']['issues']:
                                    print(f"   • {issue}")
                                
                                confirm = input("\nContinue with this table anyway? (y/n): ").strip().lower()
                                if confirm not in ['y', 'yes']:
                                    continue
                            else:
                                print(f"\n✅ Table '{selected_info['name']}' is ready for docking!")
                                print(f"   🚀 {selected_info['sdf_blob_count']:,} molecules with 3D coordinates available")
                            
                            return selected_info['name']
                        else:
                            print(f"❌ Invalid selection. Please enter 1-{len(table_data)}")
                            continue
                            
                    except ValueError:
                        # Try as table name
                        matching_tables = [t for t in table_data if t['name'].lower() == selection.lower()]
                        if matching_tables:
                            return matching_tables[0]['name']
                        else:
                            print(f"❌ Table '{selection}' not found")
                            continue
                            
                except KeyboardInterrupt:
                    print("\n❌ Selection cancelled")
                    return None
                    
        except Exception as e:
            print(f"❌ Error in table selection: {e}")
            return None


    def _show_table_details_for_docking(self, table_info: Dict[str, Any]) -> None:
        """
        Show detailed information about a table for docking preparation.
        
        Args:
            table_info (Dict): Table information dictionary
        """
        try:
            print(f"\n📋 DOCKING TABLE DETAILS: '{table_info['name']}'")
            print("=" * 60)
            print(f"🏷️  Type: {table_info['type']}")
            print(f"📊 Total Compounds: {table_info['compound_count']:,}")
            print(f"🧪 Has SMILES: {'Yes' if table_info['has_smiles'] else 'No'}")
            print(f"✅ Valid SMILES: {table_info['valid_smiles_count']} (sample)")
            
            # SDF blob information
            print(f"🔬 Has 3D Coordinates: {'Yes' if table_info['has_sdf_blob'] else 'No'}")
            if table_info['has_sdf_blob']:
                print(f"🚀 Prepared Molecules: {table_info['sdf_blob_count']:,}")
                if table_info['sdf_blob_count'] > 0:
                    prep_percentage = (table_info['sdf_blob_count'] / table_info['compound_count']) * 100
                    print(f"📈 Preparation Coverage: {prep_percentage:.1f}%")
            
            readiness = table_info['docking_readiness']
            print(f"\n🎯 DOCKING READINESS:")
            print(f"   Status: {readiness['docking_status']}")
            print(f"   Quality Score: {readiness['score']}/100")
            print(f"   Overall Status: {readiness['status']}")
            
            if readiness['issues']:
                print(f"   Issues:")
                for issue in readiness['issues']:
                    print(f"      • {issue}")
            else:
                print(f"   ✅ No issues detected")
            
            # Enhanced recommendations based on docking readiness
            print(f"\n💡 RECOMMENDATIONS:")
            if table_info['docking_status'] == 'ready for docking':
                print("   ✅ This table is ready for docking experiments")
                print("   🚀 You can proceed with receptor preparation and docking setup")
                print("   📊 All molecules have prepared 3D coordinates")
            else:
                print("   🔧 This table needs 3D molecular preparation before docking")
                
                if not table_info['has_sdf_blob']:
                    print("   📝 Steps needed:")
                    print("      1. Generate 3D conformers: chemspace.generate_mols_in_table()")
                    print("      2. This will create SDF blobs with 3D coordinates")
                    print("      3. Table will then be ready for docking")
                elif table_info['sdf_blob_count'] == 0:
                    print("   📝 SDF blob column exists but no molecules prepared")
                    print("      • Run: chemspace.generate_mols_in_table() to populate 3D coordinates")
                else:
                    partial_prep = (table_info['sdf_blob_count'] / table_info['compound_count']) * 100
                    print(f"   📊 Partial preparation: {partial_prep:.1f}% molecules ready")
                    print("      • Some molecules have 3D coordinates, others need preparation")
            
            print("=" * 60)
            
        except Exception as e:
            print(f"❌ Error showing table details: {e}")