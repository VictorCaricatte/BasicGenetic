"""
# ==============================================================================
# Author:       Victor S Caricatte De Araújo
# Email:        victorleniwys@gmail.com or victorsc@ufmg.br
# Intitution:   Universidade federal de Minas Gerais
# Version:      1.3.9
# Date:         Jun, 4
# ...................................
# ==============================================================================
"""

import sys
import os
import json
import csv
import gzip
from typing import Optional, Dict, List, Tuple, Set
from collections import defaultdict
import time
import concurrent.futures
import ftplib
from xml.etree import ElementTree
import pandas as pd
import numpy as np
import vobject
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
                            QLabel, QPushButton, QFileDialog, QComboBox, QCheckBox, 
                            QSlider, QTableWidget, QTableWidgetItem, QTabWidget, 
                            QMessageBox, QSplitter, QLineEdit, QTextEdit, QProgressBar)
from PyQt6.QtCore import Qt, QSize, QThread, pyqtSignal
from PyQt6.QtGui import QColor
import pyqtgraph as pg
from pyqtgraph import PlotWidget, BarGraphItem, HistogramLUTItem


class DatabaseWorker(QThread):
    """Thread for database operations to prevent UI freezing"""
    progress = pyqtSignal(int)
    message = pyqtSignal(str)
    finished = pyqtSignal(dict)

    def __init__(self, variants, db_type):
        super().__init__()
        self.variants = variants
        self.db_type = db_type
        self.annotator = VariantAnnotator()

    def run(self):
        results = {}
        total = len(self.variants)
        for i, variant in enumerate(self.variants):
            if self.db_type == 'all':
                results[f"{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0]}"] = \
                    self.annotator.annotate_variant(str(variant.CHROM), variant.POS, 
                                                   variant.REF, str(variant.ALT[0]))
            elif self.db_type == 'clinvar':
                results[f"{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0]}"] = \
                    self.annotator.query_clinvar(str(variant.CHROM), variant.POS, 
                                               variant.REF, str(variant.ALT[0]))
            # Similar for other db_types
            
            self.progress.emit(int((i+1)/total * 100))
            self.message.emit(f"Processing variant {i+1}/{total}")
        
        self.finished.emit(results)


class VariantAnnotator:
    """Handles variant annotation from various databases"""
    
    def __init__(self):
        self.dbsnp_cache = {}
        self.clinvar_cache = {}
        self.gnomad_cache = {}
        self.gene_cache = {}
        
        # Initialize BigQuery client if credentials exist
        self.bq_client = None
        if os.path.exists('google_credentials.json'):
            self.bq_client = bigquery.Client.from_service_account_json('google_credentials.json')
    
    def annotate_variant(self, chrom: str, pos: int, ref: str, alt: str) -> Dict:
        """Annotate a single variant with information from all databases"""
        variant_id = f"{chrom}_{pos}_{ref}_{alt}"
        
        annotation = {
            'dbsnp': self.query_dbsnp(chrom, pos, ref, alt),
            'clinvar': self.query_clinvar(chrom, pos, ref, alt),
            'gnomad': self.query_gnomad(chrom, pos, ref, alt),
            'gene_info': self.get_gene_info(chrom, pos)
        }
        
        return annotation
    
    def query_dbsnp(self, chrom: str, pos: int, ref: str, alt: str) -> Dict:
        """Query dbSNP via Entrez e-utilities"""
        variant_id = f"{chrom}_{pos}_{ref}_{alt}"
        
        if variant_id in self.dbsnp_cache:
            return self.dbsnp_cache[variant_id]
        
        # Mock implementation - in real app, use Bio.Entrez or requests
        # to query NCBI's e-utilities API
        time.sleep(0.1)  # Simulate network delay
        
        # Example response structure
        result = {
            'rsid': f"rs{np.random.randint(1000000, 9999999)}",
            'allele_frequencies': {
                '1000G': np.random.uniform(0, 0.5),
                'TOPMED': np.random.uniform(0, 0.5)
            },
            'clinical_significance': 'Uncertain significance' if np.random.random() > 0.7 else 'None',
            'gene_consequence': 'missense_variant' if np.random.random() > 0.5 else 'synonymous_variant'
        }
        
        self.dbsnp_cache[variant_id] = result
        return result
    
    def query_clinvar(self, chrom: str, pos: int, ref: str, alt: str) -> Dict:
        """Query ClinVar via FTP or API"""
        variant_id = f"{chrom}_{pos}_{ref}_{alt}"
        
        if variant_id in self.clinvar_cache:
            return self.clinvar_cache[variant_id]
        
        # Mock implementation - in real app, use FTP or API
        time.sleep(0.1)  # Simulate network delay
        
        # Example response structure
        significance = ['Pathogenic', 'Likely pathogenic', 'Uncertain significance',
                      'Likely benign', 'Benign'][np.random.randint(0, 5)]
        
        result = {
            'clinical_significance': significance,
            'review_status': f"{np.random.randint(1, 3)} stars",
            'conditions': ['Breast cancer'] if significance in ['Pathogenic', 'Likely pathogenic'] else ['None'],
            'submitter': 'ClinVar' if np.random.random() > 0.5 else 'GeneDx'
        }
        
        self.clinvar_cache[variant_id] = result
        return result
    
    def query_gnomad(self, chrom: str, pos: int, ref: str, alt: str) -> Dict:
        """Query gnomAD via BigQuery or local files"""
        variant_id = f"{chrom}_{pos}_{ref}_{alt}"
        
        if variant_id in self.gnomad_cache:
            return self.gnomad_cache[variant_id]
        
        if self.bq_client:
            # Real BigQuery implementation
            query = f"""
                SELECT 
                    AF, AF_afr, AF_amr, AF_asj, AF_eas, AF_fin, AF_nfe, AF_oth,
                    popmax, AC, AN
                FROM `gnomad-public-data.gnomad_genomes.VCF_FIELDS`
                WHERE chromosome = '{chrom}' 
                AND position = {pos}
                AND reference = '{ref}'
                AND alternate = '{alt}'
                LIMIT 1
            """
            
            try:
                query_job = self.bq_client.query(query)
                results = query_job.result()
                
                if results.total_rows > 0:
                    row = next(results)
                    result = {
                        'af': row['AF'],
                        'af_popmax': row['popmax'],
                        'ac': row['AC'],
                        'an': row['AN'],
                        'af_continental': {
                            'afr': row['AF_afr'],
                            'amr': row['AF_amr'],
                            'asj': row['AF_asj'],
                            'eas': row['AF_eas'],
                            'fin': row['AF_fin'],
                            'nfe': row['AF_nfe'],
                            'oth': row['AF_oth']
                        }
                    }
                    self.gnomad_cache[variant_id] = result
                    return result
            except Exception as e:
                print(f"Error querying gnomAD: {e}")
        
        # Fallback to mock data if BigQuery fails or not configured
        time.sleep(0.1)  # Simulate network delay
        
        result = {
            'af': np.random.uniform(0, 0.1),
            'af_popmax': np.random.uniform(0, 0.2),
            'ac': np.random.randint(1, 100),
            'an': np.random.randint(1000, 10000),
            'af_continental': {
                'afr': np.random.uniform(0, 0.1),
                'amr': np.random.uniform(0, 0.1),
                'asj': np.random.uniform(0, 0.1),
                'eas': np.random.uniform(0, 0.1),
                'fin': np.random.uniform(0, 0.1),
                'nfe': np.random.uniform(0, 0.1),
                'oth': np.random.uniform(0, 0.1)
            }
        }
        
        self.gnomad_cache[variant_id] = result
        return result
    
    def get_gene_info(self, chrom: str, pos: int) -> Dict:
        """Get gene information for a genomic position"""
        pos_key = f"{chrom}_{pos}"
        
        if pos_key in self.gene_cache:
            return self.gene_cache[pos_key]
        
        # Mock implementation - in real app, use Ensembl API or local database
        time.sleep(0.05)  # Simulate network delay
        
        genes = ['BRCA1', 'BRCA2', 'TP53', 'EGFR', 'CFTR']
        consequences = ['missense_variant', 'synonymous_variant', 'stop_gained',
                       'splice_acceptor_variant', 'intron_variant']
        
        result = {
            'gene': genes[np.random.randint(0, len(genes))],
            'transcript': f"ENST0000{np.random.randint(1000000, 9999999)}",
            'consequence': consequences[np.random.randint(0, len(consequences))],
            'impact_score': np.random.randint(1, 5),
            'protein_change': f"p.{'ACDEFGHIKLMNPQRSTVWY'[np.random.randint(0, 20)]}"
                              f"{np.random.randint(1, 500)}"
                              f"{'ACDEFGHIKLMNPQRSTVWY'[np.random.randint(0, 20)]}"
        }
        
        self.gene_cache[pos_key] = result
        return result


class VCFLoader:
    """Handles loading and parsing VCF files"""
    
    def __init__(self):
        self.vcf_reader = None
        self.variants = []
        self.samples = []
        self.file_path = ""
        self.file_name = ""
    
    def load_vcf(self, file_path: str) -> bool:
        """Load a VCF file and parse its contents"""
        try:
            self.file_path = file_path
            self.file_name = os.path.basename(file_path)
            self.vcf_reader = vcf.Reader(filename=file_path)
            self.variants = list(self.vcf_reader)
            self.samples = self.vcf_reader.samples
            return True
        except Exception as e:
            print(f"Error loading VCF file: {e}")
            return False
    
    def get_variant_count(self) -> int:
        """Return the number of variants loaded"""
        return len(self.variants)
    
    def get_sample_names(self) -> List[str]:
        """Return the list of sample names"""
        return self.samples
    
    def get_filtered_vcf(self, variant_indices: List[int]) -> str:
        """Generate a reduced VCF containing only the specified variants"""
        if not self.vcf_reader:
            return ""
        
        output = []
        
        # Write header
        for header_line in self.vcf_reader._header_lines:
            output.append(header_line)
        
        # Write selected variants
        for idx in variant_indices:
            if 0 <= idx < len(self.variants):
                variant = self.variants[idx]
                output.append("\t".join(map(str, [
                    variant.CHROM, variant.POS, variant.ID or '.',
                    variant.REF, ",".join(map(str, variant.ALT)),
                    variant.QUAL or '.', ";".join(variant.FILTER) if variant.FILTER else 'PASS',
                    self._format_info(variant), self._format_samples(variant)
                ])))
        
        return "\n".join(output)
    
    def _format_info(self, variant) -> str:
        """Format INFO field for VCF output"""
        if not hasattr(variant, 'INFO'):
            return '.'
        
        info_parts = []
        for key, value in variant.INFO.items():
            if isinstance(value, bool):
                info_parts.append(key)
            elif isinstance(value, (list, tuple)):
                info_parts.append(f"{key}={','.join(map(str, value))}")
            else:
                info_parts.append(f"{key}={value}")
        
        return ";".join(info_parts)
    
    def _format_samples(self, variant) -> str:
        """Format SAMPLE fields for VCF output"""
        if not self.samples:
            return '.'
        
        format_fields = variant.FORMAT.split(':') if hasattr(variant, 'FORMAT') else []
        sample_data = []
        
        for sample in self.samples:
            sample_info = variant.genotype(sample)
            sample_parts = []
            
            for field in format_fields:
                if hasattr(sample_info.data, field):
                    value = getattr(sample_info.data, field)
                    if isinstance(value, (list, tuple)):
                        sample_parts.append(",".join(map(str, value)))
                    else:
                        sample_parts.append(str(value))
                else:
                    sample_parts.append('.')
            
            sample_data.append(":".join(sample_parts))
        
        return "\t".join(sample_data)


class MultiVCFLoader:
    """Handles loading and comparing multiple VCF files"""
    
    def __init__(self):
        self.vcf_loaders = []
        self.combined_variants = []
        self.variant_sources = {}  # Maps variant keys to source file indices
    
    def add_vcf(self, file_path: str) -> bool:
        """Add a VCF file to the collection"""
        loader = VCFLoader()
        if loader.load_vcf(file_path):
            self.vcf_loaders.append(loader)
            return True
        return False
    
    def compare_vcfs(self) -> Dict[str, Set]:
        """Compare the loaded VCFs and find intersections/differences"""
        variant_sets = []
        
        # Create sets of variant keys for each VCF
        for loader in self.vcf_loaders:
            variant_keys = set()
            for variant in loader.variants:
                key = self._variant_key(variant)
                variant_keys.add(key)
            variant_sets.append(variant_keys)
        
        if not variant_sets:
            return {}
        
        # Find common and unique variants
        comparisons = {
            'common_to_all': set(variant_sets[0]),
            'unique_to_each': [set() for _ in range(len(variant_sets))]
        }
        
        for i, variant_set in enumerate(variant_sets):
            comparisons['common_to_all'] &= variant_set
            comparisons['unique_to_each'][i] = variant_set - comparisons['common_to_all']
            
            for j, other_set in enumerate(variant_sets):
                if i != j:
                    comparisons['unique_to_each'][i] -= other_set
        
        return comparisons
    
    def _variant_key(self, variant) -> str:
        """Create a unique key for a variant"""
        return f"{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0]}"
    
    def get_loader(self, index: int) -> Optional[VCFLoader]:
        """Get a VCFLoader by index"""
        if 0 <= index < len(self.vcf_loaders):
            return self.vcf_loaders[index]
        return None


class VariantTableModel:
    """Manages variant data for display in tables and plots"""
    
    def __init__(self):
        self.variants = []
        self.filtered_variants = []
        self.current_filters = {}
        self.variant_annotations = {}  # Cache for variant annotations
        self.vcf_loader = None
    
    def set_variants(self, variants: List, vcf_loader: Optional[VCFLoader] = None):
        """Set the variants to be managed"""
        self.variants = variants
        self.filtered_variants = variants.copy()
        self.current_filters = {}
        self.vcf_loader = vcf_loader
    
    def apply_filters(self, filters: Dict):
        """Apply filters to the variant data"""
        self.current_filters = filters
        self.filtered_variants = self.variants.copy()
        
        if 'impact' in filters and filters['impact']:
            self.filtered_variants = [
                v for v in self.filtered_variants 
                if self._get_impact_type(v) in filters['impact']
            ]
        
        if 'max_freq' in filters and filters['max_freq'] is not None:
            self.filtered_variants = [
                v for v in self.filtered_variants 
                if self._get_variant_frequency(v) <= filters['max_freq']
            ]
        
        if 'chromosome' in filters and filters['chromosome']:
            self.filtered_variants = [
                v for v in self.filtered_variants 
                if str(v.CHROM) == filters['chromosome']
            ]
        
        if 'gene_search' in filters and filters['gene_search']:
            search_term = filters['gene_search'].lower()
            self.filtered_variants = [
                v for v in self.filtered_variants
                if self._gene_matches(v, search_term)
            ]
        
        if 'variant_search' in filters and filters['variant_search']:
            search_term = filters['variant_search'].lower()
            self.filtered_variants = [
                v for v in self.filtered_variants
                if self._variant_matches(v, search_term)
            ]
        
        if 'region_search' in filters and filters['region_search']:
            chrom, start, end = self._parse_region(filters['region_search'])
            if chrom:
                self.filtered_variants = [
                    v for v in self.filtered_variants
                    if str(v.CHROM) == chrom and start <= v.POS <= end
                ]
    
    def _get_impact_type(self, variant) -> str:
        """Determine the impact type of a variant"""
        if not hasattr(variant, 'INFO') or 'ANN' not in variant.INFO:
            return 'unknown'
        
        ann = variant.INFO['ANN'][0]
        if 'missense' in ann:
            return 'missense'
        elif 'nonsense' in ann:
            return 'nonsense'
        elif 'synonymous' in ann:
            return 'synonymous'
        return 'other'
    
    def _get_variant_frequency(self, variant) -> float:
        """Get the population frequency of a variant"""
        if hasattr(variant, 'INFO') and 'AF' in variant.INFO:
            return float(variant.INFO['AF'][0])
        return 0.0
    
    def _gene_matches(self, variant, search_term: str) -> bool:
        """Check if variant matches gene search"""
        # In a real implementation, this would check gene annotations
        return search_term in str(variant.CHROM).lower() or \
               search_term in str(variant.POS).lower()  # Simplified
    
    def _variant_matches(self, variant, search_term: str) -> bool:
        """Check if variant matches ID search"""
        return (variant.ID and search_term in variant.ID.lower()) or \
               search_term in f"{variant.CHROM}_{variant.POS}".lower()
    
    def _parse_region(self, region_str: str) -> Tuple[Optional[str], int, int]:
        """Parse a genomic region string (e.g., 'chr1:1000-2000')"""
        if ':' not in region_str or '-' not in region_str:
            return None, 0, 0
        
        try:
            chrom_part, pos_part = region_str.split(':')
            start_str, end_str = pos_part.split('-')
            return chrom_part, int(start_str), int(end_str)
        except (ValueError, AttributeError):
            return None, 0, 0
    
    def get_filtered_variants(self) -> List:
        """Get the currently filtered variants"""
        return self.filtered_variants
    
    def get_filtered_indices(self) -> List[int]:
        """Get the indices of filtered variants in the original list"""
        variant_set = {id(v): idx for idx, v in enumerate(self.variants)}
        return [variant_set[id(v)] for v in self.filtered_variants]
    
    def get_chromosome_counts(self) -> Dict[str, int]:
        """Get counts of variants per chromosome"""
        chrom_counts = {}
        for v in self.variants:
            chrom = str(v.CHROM)
            chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        return chrom_counts
    
    def get_impact_counts(self) -> Dict[str, int]:
        """Get counts of variants by impact type"""
        impact_counts = {'missense': 0, 'nonsense': 0, 'synonymous': 0, 'other': 0, 'unknown': 0}
        for v in self.variants:
            impact = self._get_impact_type(v)
            impact_counts[impact] += 1
        return impact_counts
    
    def get_quality_distribution(self) -> np.ndarray:
        """Get quality scores for histogram"""
        quals = []
        for v in self.variants:
            if v.QUAL is not None:
                quals.append(v.QUAL)
        return np.array(quals)
    
    def get_frequency_distribution(self) -> np.ndarray:
        """Get frequency values for histogram"""
        freqs = []
        for v in self.variants:
            if hasattr(v, 'INFO') and 'AF' in v.INFO:
                freqs.append(float(v.INFO['AF'][0]))
        return np.array(freqs)
    
    def save_to_file(self, file_path: str, format_type: str = 'vcf') -> bool:
        """Save filtered variants to a file"""
        if not self.vcf_loader:
            return False
        
        try:
            if format_type == 'vcf':
                with open(file_path, 'w') as f:
                    f.write(self.vcf_loader.get_filtered_vcf(self.get_filtered_indices()))
            elif format_type == 'csv':
                with open(file_path, 'w', newline='') as f:
                    writer = csv.writer(f)
                    # Write header
                    writer.writerow([
                        'Chromosome', 'Position', 'ID', 'Reference', 'Alternate',
                        'Quality', 'Filter', 'Impact', 'Gene', 'Consequence'
                    ])
                    
                    # Write data
                    for variant in self.filtered_variants:
                        gene_info = self.variant_annotations.get(
                            f"{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0]}", {}
                        ).get('gene_info', {})
                        
                        writer.writerow([
                            str(variant.CHROM), variant.POS, variant.ID or '',
                            variant.REF, ",".join(map(str, variant.ALT)),
                            variant.QUAL or '', ";".join(variant.FILTER) if variant.FILTER else 'PASS',
                            self._get_impact_type(variant),
                            gene_info.get('gene', ''),
                            gene_info.get('consequence', '')
                        ])
            elif format_type == 'json':
                data = []
                for variant in self.filtered_variants:
                    var_data = {
                        'chromosome': str(variant.CHROM),
                        'position': variant.POS,
                        'id': variant.ID,
                        'reference': variant.REF,
                        'alternate': [str(a) for a in variant.ALT],
                        'quality': variant.QUAL,
                        'filter': variant.FILTER,
                        'impact': self._get_impact_type(variant),
                        'annotations': self.variant_annotations.get(
                            f"{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0]}", {}
                        )
                    }
                    data.append(var_data)
                
                with open(file_path, 'w') as f:
                    json.dump(data, f, indent=2)
            
            return True
        except Exception as e:
            print(f"Error saving file: {e}")
            return False


class ChromosomePlotWidget(QWidget):
    """Widget for displaying chromosome variant distribution"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('w')
        self.plot_widget.setLabel('left', 'Variant Count')
        self.plot_widget.setLabel('bottom', 'Chromosome')
        self.plot_widget.showGrid(x=True, y=True)
        
        self.bar_graph = None
        
        layout = QVBoxLayout()
        layout.addWidget(self.plot_widget)
        self.setLayout(layout)
    
    def update_plot(self, chrom_counts: Dict[str, int]):
        """Update the plot with new chromosome data"""
        chromosomes = sorted(chrom_counts.keys(), key=lambda x: int(x) if x.isdigit() else float('inf'))
        counts = [chrom_counts[chrom] for chrom in chromosomes]
        
        x = np.arange(len(chromosomes))
        y = np.array(counts)
        
        if self.bar_graph is not None:
            self.plot_widget.removeItem(self.bar_graph)
        
        # Create bars with different colors
        self.bar_graph = BarGraphItem(x=x, height=y, width=0.6, brush=pg.mkBrush('#3498db'))
        self.plot_widget.addItem(self.bar_graph)
        
        # Set x-axis ticks
        self.plot_widget.getPlotItem().getAxis('bottom').setTicks([[(i, chrom) for i, chrom in enumerate(chromosomes)]])


class ImpactPlotWidget(QWidget):
    """Widget for displaying variant impact distribution"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('w')
        self.plot_widget.setLabel('left', 'Variant Count')
        self.plot_widget.setLabel('bottom', 'Impact Type')
        self.plot_widget.showGrid(x=True, y=True)
        
        self.bar_graph = None
        
        layout = QVBoxLayout()
        layout.addWidget(self.plot_widget)
        self.setLayout(layout)
    
    def update_plot(self, impact_counts: Dict[str, int]):
        """Update the plot with new impact data"""
        impacts = ['missense', 'nonsense', 'synonymous', 'other', 'unknown']
        counts = [impact_counts.get(imp, 0) for imp in impacts]
        
        x = np.arange(len(impacts))
        y = np.array(counts)
        
        if self.bar_graph is not None:
            self.plot_widget.removeItem(self.bar_graph)
        
        # Create colored bars
        colors = ['#e74c3c', '#f39c12', '#2ecc71', '#9b59b6', '#34495e']
        self.bar_graph = BarGraphItem(x=x, height=y, width=0.6, 
                                    brushes=[pg.mkBrush(color) for color in colors])
        self.plot_widget.addItem(self.bar_graph)
        
        # Set x-axis ticks
        self.plot_widget.getPlotItem().getAxis('bottom').setTicks([[(i, imp) for i, imp in enumerate(impacts)]])


class QualityHistogramWidget(QWidget):
    """Widget for displaying quality score distribution"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('w')
        self.plot_widget.setLabel('left', 'Count')
        self.plot_widget.setLabel('bottom', 'Quality Score')
        self.plot_widget.showGrid(x=True, y=True)
        
        self.histogram = None
        
        layout = QVBoxLayout()
        layout.addWidget(self.plot_widget)
        self.setLayout(layout)
    
    def update_plot(self, quality_scores: np.ndarray):
        """Update the plot with new quality data"""
        if self.histogram is not None:
            self.plot_widget.removeItem(self.histogram)
        
        if len(quality_scores) == 0:
            return
        
        y, x = np.histogram(quality_scores, bins=50)
        self.histogram = pg.BarGraphItem(x=x[:-1], height=y, width=0.8*(x[1]-x[0]), brush=pg.mkBrush('#3498db'))
        self.plot_widget.addItem(self.histogram)
        
        # Add a line at Q30 (common quality threshold)
        q30_line = pg.InfiniteLine(pos=30, angle=90, pen=pg.mkPen('#e74c3c', width=2))
        self.plot_widget.addItem(q30_line)
        
        # Add text label for Q30
        q30_text = pg.TextItem(text="Q30", color='#e74c3c', anchor=(0, 1))
        q30_text.setPos(30, max(y)*0.9)
        self.plot_widget.addItem(q30_text)


class FrequencyHistogramWidget(QWidget):
    """Widget for displaying allele frequency distribution"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('w')
        self.plot_widget.setLabel('left', 'Count')
        self.plot_widget.setLabel('bottom', 'Allele Frequency')
        self.plot_widget.showGrid(x=True, y=True)
        
        self.histogram = None
        
        layout = QVBoxLayout()
        layout.addWidget(self.plot_widget)
        self.setLayout(layout)
    
    def update_plot(self, frequencies: np.ndarray):
        """Update the plot with new frequency data"""
        if self.histogram is not None:
            self.plot_widget.removeItem(self.histogram)
        
        if len(frequencies) == 0:
            return
        
        y, x = np.histogram(frequencies, bins=50)
        self.histogram = pg.BarGraphItem(x=x[:-1], height=y, width=0.8*(x[1]-x[0]), brush=pg.mkBrush('#2ecc71'))
        self.plot_widget.addItem(self.histogram)
        
        # Add a line at 1% frequency
        one_percent_line = pg.InfiniteLine(pos=0.01, angle=90, pen=pg.mkPen('#e74c3c', width=2))
        self.plot_widget.addItem(one_percent_line)
        
        # Add text label for 1%
        one_percent_text = pg.TextItem(text="1%", color='#e74c3c', anchor=(0, 1))
        one_percent_text.setPos(0.01, max(y)*0.9)
        self.plot_widget.addItem(one_percent_text)


class ComparisonPlotWidget(QWidget):
    """Widget for displaying comparison between multiple VCFs"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('w')
        self.plot_widget.setLabel('left', 'Variant Count')
        self.plot_widget.setLabel('bottom', 'Sample')
        self.plot_widget.showGrid(x=True, y=True)
        
        self.bar_graphs = []
        
        layout = QVBoxLayout()
        layout.addWidget(self.plot_widget)
        self.setLayout(layout)
    
    def update_plot(self, comparisons: Dict[str, List], sample_names: List[str]):
        """Update the plot with comparison data"""
        # Clear previous bars
        for bar in self.bar_graphs:
            self.plot_widget.removeItem(bar)
        self.bar_graphs = []
        
        if not comparisons or not sample_names:
            return
        
        # Prepare data for stacked bars
        unique_counts = [len(s) for s in comparisons['unique_to_each']]
        common_count = len(comparisons['common_to_all'])
        
        x = np.arange(len(sample_names))
        y_unique = np.array(unique_counts)
        y_common = np.array([common_count] * len(sample_names))
        
        # Create stacked bars
        common_bars = BarGraphItem(x=x, height=y_common, width=0.6, brush=pg.mkBrush('#3498db'))
        unique_bars = BarGraphItem(x=x, height=y_unique, width=0.6, brush=pg.mkBrush('#e74c3c'))
        
        self.plot_widget.addItem(common_bars)
        self.plot_widget.addItem(unique_bars)
        self.bar_graphs.extend([common_bars, unique_bars])
        
        # Add legend
        legend = pg.LegendItem(offset=(10, 10))
        legend.addItem(common_bars, 'Common to all')
        legend.addItem(unique_bars, 'Unique to sample')
        legend.setParentItem(self.plot_widget.getPlotItem())
        self.bar_graphs.append(legend)
        
        # Set x-axis ticks
        self.plot_widget.getPlotItem().getAxis('bottom').setTicks([[(i, name) for i, name in enumerate(sample_names)]])


class VariantTableWidget(QTableWidget):
    """Table widget for displaying variant information"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setColumnCount(10)
        self.setHorizontalHeaderLabels([
            'Chromosome', 'Position', 'ID', 'Ref', 'Alt', 
            'Quality', 'Filter', 'Impact', 'Gene', 'Consequence'
        ])
        self.setSortingEnabled(True)
        self.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        
        # Adjust column widths
        self.setColumnWidth(0, 100)  # Chromosome
        self.setColumnWidth(1, 80)   # Position
        self.setColumnWidth(2, 120)  # ID
        self.setColumnWidth(3, 60)   # Ref
        self.setColumnWidth(4, 60)   # Alt
        self.setColumnWidth(5, 80)   # Quality
        self.setColumnWidth(6, 100)  # Filter
        self.setColumnWidth(7, 100)  # Impact
        self.setColumnWidth(8, 100)  # Gene
        self.setColumnWidth(9, 150)  # Consequence
    
    def update_table(self, variants: List, annotations: Optional[Dict] = None):
        """Update the table with new variant data"""
        self.setRowCount(len(variants))
        
        for row, variant in enumerate(variants):
            variant_key = f"{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0]}"
            gene_info = annotations.get(variant_key, {}).get('gene_info', {}) if annotations else {}
            
            # Chromosome
            self.setItem(row, 0, QTableWidgetItem(str(variant.CHROM)))
            
            # Position
            self.setItem(row, 1, QTableWidgetItem(str(variant.POS)))
            
            # ID
            id_item = QTableWidgetItem(variant.ID if variant.ID else '')
            self.setItem(row, 2, id_item)
            
            # Reference allele
            self.setItem(row, 3, QTableWidgetItem(variant.REF))
            
            # Alternate allele(s)
            alt_text = ','.join(str(a) for a in variant.ALT)
            self.setItem(row, 4, QTableWidgetItem(alt_text))
            
            # Quality
            qual_item = QTableWidgetItem(str(variant.QUAL) if variant.QUAL else '')
            self.setItem(row, 5, qual_item)
            
            # Filter
            filter_text = ';'.join(variant.FILTER) if variant.FILTER else 'PASS'
            filter_item = QTableWidgetItem(filter_text)
            
            # Color filter status
            if filter_text != 'PASS':
                filter_item.setBackground(QColor(255, 200, 200))
            self.setItem(row, 6, filter_item)
            
            # Impact (simplified)
            impact = self._get_impact_type(variant)
            impact_item = QTableWidgetItem(impact)
            
            # Color by impact type
            if impact == 'missense':
                impact_item.setBackground(QColor(255, 200, 200))
            elif impact == 'nonsense':
                impact_item.setBackground(QColor(255, 150, 150))
            elif impact == 'synonymous':
                impact_item.setBackground(QColor(200, 255, 200))
                
            self.setItem(row, 7, impact_item)
            
            # Gene
            gene_item = QTableWidgetItem(gene_info.get('gene', ''))
            self.setItem(row, 8, gene_item)
            
            # Consequence
            consequence_item = QTableWidgetItem(gene_info.get('consequence', ''))
            self.setItem(row, 9, consequence_item)
    
    def _get_impact_type(self, variant) -> str:
        """Determine the impact type of a variant (simplified)"""
        if not hasattr(variant, 'INFO') or 'ANN' not in variant.INFO:
            return 'unknown'
        
        ann = variant.INFO['ANN'][0]
        if 'missense' in ann:
            return 'missense'
        elif 'nonsense' in ann:
            return 'nonsense'
        elif 'synonymous' in ann:
            return 'synonymous'
        return 'other'


class VariantDetailWidget(QWidget):
    """Widget for displaying detailed information about a selected variant"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumWidth(300)
        
        layout = QVBoxLayout()
        
        # Title
        self.title_label = QLabel("Variant Details")
        self.title_label.setStyleSheet("font-weight: bold; font-size: 14px;")
        layout.addWidget(self.title_label)
        
        # Details table
        self.detail_table = QTableWidget()
        self.detail_table.setColumnCount(2)
        self.detail_table.setHorizontalHeaderLabels(['Field', 'Value'])
        self.detail_table.verticalHeader().setVisible(False)
        self.detail_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.detail_table.horizontalHeader().setStretchLastSection(True)
        self.detail_table.setColumnWidth(0, 120)
        layout.addWidget(self.detail_table)
        
        # Gene info section
        self.gene_label = QLabel("Gene Information")
        self.gene_label.setStyleSheet("font-weight: bold; font-size: 14px; margin-top: 10px;")
        layout.addWidget(self.gene_label)
        
        self.gene_table = QTableWidget()
        self.gene_table.setColumnCount(2)
        self.gene_table.setHorizontalHeaderLabels(['Field', 'Value'])
        self.gene_table.verticalHeader().setVisible(False)
        self.gene_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.gene_table.horizontalHeader().setStretchLastSection(True)
        self.gene_table.setColumnWidth(0, 120)
        layout.addWidget(self.gene_table)
        
        # Annotation section
        self.annotation_label = QLabel("Database Annotations")
        self.annotation_label.setStyleSheet("font-weight: bold; font-size: 14px; margin-top: 10px;")
        layout.addWidget(self.annotation_label)
        
        self.annotation_table = QTableWidget()
        self.annotation_table.setColumnCount(2)
        self.annotation_table.setHorizontalHeaderLabels(['Database', 'Annotation'])
        self.annotation_table.verticalHeader().setVisible(False)
        self.annotation_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.annotation_table.horizontalHeader().setStretchLastSection(True)
        self.annotation_table.setColumnWidth(0, 120)
        layout.addWidget(self.annotation_table)
        
        layout.addStretch()
        self.setLayout(layout)
    
    def update_details(self, variant, annotation: Optional[Dict] = None):
        """Update the details with information about the selected variant"""
        self.detail_table.setRowCount(0)
        self.gene_table.setRowCount(0)
        self.annotation_table.setRowCount(0)
        
        if variant is None:
            return
        
        # Basic variant information
        basic_info = [
            ('Chromosome', str(variant.CHROM)),
            ('Position', str(variant.POS)),
            ('ID', variant.ID if variant.ID else 'None'),
            ('Reference', variant.REF),
            ('Alternate', ', '.join(str(a) for a in variant.ALT)),
            ('Quality', str(variant.QUAL) if variant.QUAL else 'None'),
            ('Filter', '; '.join(variant.FILTER) if variant.FILTER else 'PASS'),
            ('Sample Count', str(len(variant.samples))),
        ]
        
        self.detail_table.setRowCount(len(basic_info))
        for row, (field, value) in enumerate(basic_info):
            self.detail_table.setItem(row, 0, QTableWidgetItem(field))
            self.detail_table.setItem(row, 1, QTableWidgetItem(value))
        
        # Add INFO fields if present
        if hasattr(variant, 'INFO') and variant.INFO:
            info_start_row = self.detail_table.rowCount()
            self.detail_table.setRowCount(info_start_row + len(variant.INFO))
            
            for row, (key, value) in enumerate(variant.INFO.items(), start=info_start_row):
                self.detail_table.setItem(row, 0, QTableWidgetItem(f"INFO.{key}"))
                
                # Handle different types of INFO values
                if isinstance(value, list):
                    val_str = ', '.join(str(v) for v in value)
                else:
                    val_str = str(value)
                
                self.detail_table.setItem(row, 1, QTableWidgetItem(val_str))
        
        # Update gene information
        if annotation and 'gene_info' in annotation:
            gene_info = annotation['gene_info']
            gene_fields = [
                ('Gene', gene_info.get('gene', '')),
                ('Transcript', gene_info.get('transcript', '')),
                ('Consequence', gene_info.get('consequence', '')),
                ('Impact Score', str(gene_info.get('impact_score', ''))),
                ('Protein Change', gene_info.get('protein_change', ''))
            ]
            
            self.gene_table.setRowCount(len(gene_fields))
            for row, (field, value) in enumerate(gene_fields):
                self.gene_table.setItem(row, 0, QTableWidgetItem(field))
                self.gene_table.setItem(row, 1, QTableWidgetItem(value))
        
        # Update annotation information
        if annotation:
            db_annotations = {k: v for k, v in annotation.items() if k != 'gene_info'}
            self.annotation_table.setRowCount(len(db_annotations))
            
            for row, (db, data) in enumerate(db_annotations.items()):
                self.annotation_table.setItem(row, 0, QTableWidgetItem(db.upper()))
                
                # Format annotation data
                if data:
                    anno_text = '\n'.join(f"{k}: {v}" for k, v in data.items())
                else:
                    anno_text = "No annotation found"
                
                self.annotation_table.setItem(row, 1, QTableWidgetItem(anno_text))


class HelpTab(QWidget):
    """Widget for displaying help information"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        
        layout = QVBoxLayout()
        
        # Title
        title = QLabel("VariantGenV - Help")
        title.setStyleSheet("font-weight: bold; font-size: 16px;")
        layout.addWidget(title)
        
        # Introduction
        intro = QLabel("""
        This application helps visualize and analyze genetic variants from VCF files.
        Below are instructions for using the various features.
        Author: Victor Silveira Caricatte
        Universidade Federal De Minas Gerais - UFMG
        """)
        intro.setWordWrap(True)
        layout.addWidget(intro)
        
        # Sections
        sections = QTabWidget()
        
        # File Loading
        file_tab = QWidget()
        file_layout = QVBoxLayout(file_tab)
        file_layout.addWidget(QLabel("""
        <h3>Loading Files</h3>
        <p>Use <b>File > Open VCF</b> to load a single VCF file or <b>File > Add VCF</b> to load multiple files for comparison.</p>
        <p>Supported formats: uncompressed VCF (.vcf) or bgzip-compressed VCF (.vcf.gz).</p>
        <p>For large files, the application will show a progress indicator during loading.</p>
        """))
        file_layout.addStretch()
        sections.addTab(file_tab, "File Loading")
        
        # Filtering
        filter_tab = QWidget()
        filter_layout = QVBoxLayout(filter_tab)
        filter_layout.addWidget(QLabel("""
        <h3>Filtering Variants</h3>
        <p>The right panel provides several filtering options:</p>
        <ul>
            <li><b>Chromosome:</b> Filter by chromosome</li>
            <li><b>Impact Type:</b> Select which variant consequences to include</li>
            <li><b>Max Frequency:</b> Filter by allele frequency (0-100%)</li>
            <li><b>Search:</b> Find variants by gene name, variant ID, or genomic region</li>
        </ul>
        <p>Click <b>Apply Filters</b> to update the table and plots.</p>
        """))
        filter_layout.addStretch()
        sections.addTab(filter_tab, "Filtering")
        
        # Visualization
        viz_tab = QWidget()
        viz_layout = QVBoxLayout(viz_tab)
        viz_layout.addWidget(QLabel("""
        <h3>Visualization</h3>
        <p>The application provides several visualization tabs:</p>
        <ul>
            <li><b>Summary:</b> Shows distribution of variants by chromosome and impact type</li>
            <li><b>Quality/Frequency:</b> Histograms of variant quality scores and allele frequencies</li>
            <li><b>Comparison:</b> When multiple files are loaded, shows shared and unique variants</li>
            <li><b>Variants:</b> Table view of all variants with detailed information</li>
        </ul>
        <p>Click on a variant in the table to see detailed annotations in the right panel.</p>
        """))
        viz_layout.addStretch()
        sections.addTab(viz_tab, "Visualization")
        
        # Export
        export_tab = QWidget()
        export_layout = QVBoxLayout(export_tab)
        export_layout.addWidget(QLabel("""
        <h3>Exporting Data</h3>
        <p>Use <b>File > Export Variants</b> to save filtered variants in different formats:</p>
        <ul>
            <li><b>VCF:</b> Standard VCF format with only the filtered variants</li>
            <li><b>CSV:</b> Comma-separated values for spreadsheet applications</li>
            <li><b>JSON:</b> Structured format for programmatic use</li>
        </ul>
        """))
        export_layout.addStretch()
        sections.addTab(export_tab, "Exporting")
        
        layout.addWidget(sections)
        self.setLayout(layout)


class MainWindow(QMainWindow):
    """Main application window"""
    
    def __init__(self):
        super().__init__()
        
        # Initialize components
        self.vcf_loader = VCFLoader()
        self.multi_loader = MultiVCFLoader()
        self.variant_annotator = VariantAnnotator()
        self.variant_model = VariantTableModel()
        self.current_file_index = 0
        
        # Set up UI
        self.setWindowTitle("VariantGenV")
        self.setGeometry(100, 100, 1200, 800)
        
        # Create central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        main_layout = QHBoxLayout(central_widget)
        
        # Create splitter for main content
        splitter = QSplitter(Qt.Orientation.Horizontal)
        
        # Left panel (plots and table)
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        
        # Create tab widget for different views
        self.tab_widget = QTabWidget()
        
        # Summary tab
        summary_tab = QWidget()
        summary_layout = QVBoxLayout(summary_tab)
        
        # Chromosome plot
        self.chrom_plot = ChromosomePlotWidget()
        summary_layout.addWidget(QLabel("Variants by Chromosome"))
        summary_layout.addWidget(self.chrom_plot)
        
        # Impact plot
        self.impact_plot = ImpactPlotWidget()
        summary_layout.addWidget(QLabel("Variants by Impact Type"))
        summary_layout.addWidget(self.impact_plot)
        
        # Add summary tab
        self.tab_widget.addTab(summary_tab, "Summary")
        
        # Quality/Frequency tab
        qf_tab = QWidget()
        qf_layout = QVBoxLayout(qf_tab)
        
        # Quality histogram
        self.quality_hist = QualityHistogramWidget()
        qf_layout.addWidget(QLabel("Quality Score Distribution"))
        qf_layout.addWidget(self.quality_hist)
        
        # Frequency histogram
        self.freq_hist = FrequencyHistogramWidget()
        qf_layout.addWidget(QLabel("Allele Frequency Distribution"))
        qf_layout.addWidget(self.freq_hist)
        
        self.tab_widget.addTab(qf_tab, "Quality/Frequency")
        
        # Comparison tab
        self.comparison_tab = QWidget()
        comparison_layout = QVBoxLayout(self.comparison_tab)
        
        self.comparison_plot = ComparisonPlotWidget()
        comparison_layout.addWidget(QLabel("Variant Comparison Between Samples"))
        comparison_layout.addWidget(self.comparison_plot)
        
        self.tab_widget.addTab(self.comparison_tab, "Comparison")
        
        # Variant table tab
        table_tab = QWidget()
        table_layout = QVBoxLayout(table_tab)
        
        self.variant_table = VariantTableWidget()
        table_layout.addWidget(self.variant_table)
        
        self.tab_widget.addTab(table_tab, "Variants")
        
        # Help tab
        help_tab = HelpTab()
        self.tab_widget.addTab(help_tab, "Help")
        
        left_layout.addWidget(self.tab_widget)
        
        # Right panel (filters and details)
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        
        # Search controls
        search_group = QWidget()
        search_layout = QVBoxLayout(search_group)
        search_layout.addWidget(QLabel("Search Variants"))
        
        # Gene search
        self.gene_search = QLineEdit()
        self.gene_search.setPlaceholderText("Search by gene name...")
        search_layout.addWidget(self.gene_search)
        
        # Variant ID search
        self.variant_search = QLineEdit()
        self.variant_search.setPlaceholderText("Search by variant ID (rsID)...")
        search_layout.addWidget(self.variant_search)
        
        # Region search
        self.region_search = QLineEdit()
        self.region_search.setPlaceholderText("Search by region (e.g., chr1:10000-20000)...")
        search_layout.addWidget(self.region_search)
        
        right_layout.addWidget(search_group)
        
        # Filter controls
        filter_group = QWidget()
        filter_layout = QVBoxLayout(filter_group)
        filter_layout.addWidget(QLabel("Filters"))
        
        # Chromosome filter
        self.chrom_filter = QComboBox()
        self.chrom_filter.addItem("All Chromosomes", None)
        filter_layout.addWidget(QLabel("Chromosome:"))
        filter_layout.addWidget(self.chrom_filter)
        
        # Impact filter
        self.impact_filter_group = QWidget()
        impact_filter_layout = QVBoxLayout(self.impact_filter_group)
        impact_filter_layout.addWidget(QLabel("Impact Type:"))
        
        self.missense_check = QCheckBox("Missense")
        self.missense_check.setChecked(True)
        impact_filter_layout.addWidget(self.missense_check)
        
        self.nonsense_check = QCheckBox("Nonsense")
        self.nonsense_check.setChecked(True)
        impact_filter_layout.addWidget(self.nonsense_check)
        
        self.synonymous_check = QCheckBox("Synonymous")
        self.synonymous_check.setChecked(True)
        impact_filter_layout.addWidget(self.synonymous_check)
        
        self.other_check = QCheckBox("Other")
        self.other_check.setChecked(True)
        impact_filter_layout.addWidget(self.other_check)
        
        filter_layout.addWidget(self.impact_filter_group)
        
        # Frequency filter
        self.freq_slider = QSlider(Qt.Orientation.Horizontal)
        self.freq_slider.setRange(0, 100)
        self.freq_slider.setValue(100)
        self.freq_label = QLabel("Max Frequency: 100%")
        filter_layout.addWidget(self.freq_label)
        filter_layout.addWidget(self.freq_slider)
        
        # Apply filters button
        self.apply_filters_btn = QPushButton("Apply Filters")
        filter_layout.addWidget(self.apply_filters_btn)
        
        right_layout.addWidget(filter_group)
        
        # Variant details
        self.variant_detail = VariantDetailWidget()
        right_layout.addWidget(self.variant_detail)
        
        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        right_layout.addWidget(self.progress_bar)
        
        # Add panels to splitter
        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)
        splitter.setSizes([800, 400])
        
        main_layout.addWidget(splitter)
        
        # Menu bar
        self.create_menu_bar()
        
        # Connect signals
        self.setup_connections()
        
        # Status bar
        self.status_bar = self.statusBar()
        
        # Hide comparison tab initially
        self.tab_widget.setTabVisible(2, False)
    
    def create_menu_bar(self):
        """Create the menu bar"""
        menu_bar = self.menuBar()
        
        # File menu
        file_menu = menu_bar.addMenu("File")
        
        open_action = file_menu.addAction("Open VCF")
        open_action.triggered.connect(self.open_vcf_file)
        
        add_action = file_menu.addAction("Add VCF")
        add_action.triggered.connect(self.add_vcf_file)
        
        export_menu = file_menu.addMenu("Export Variants")
        
        export_vcf_action = export_menu.addAction("Export as VCF")
        export_vcf_action.triggered.connect(lambda: self.export_variants('vcf'))
        
        export_csv_action = export_menu.addAction("Export as CSV")
        export_csv_action.triggered.connect(lambda: self.export_variants('csv'))
        
        export_json_action = export_menu.addAction("Export as JSON")
        export_json_action.triggered.connect(lambda: self.export_variants('json'))
        
        exit_action = file_menu.addAction("Exit")
        exit_action.triggered.connect(self.close)
    
    def setup_connections(self):
        """Set up signal-slot connections"""
        # File menu
        self.apply_filters_btn.clicked.connect(self.apply_filters)
        
        # Variant table selection
        self.variant_table.itemSelectionChanged.connect(self.on_variant_selected)
        
        # Frequency slider
        self.freq_slider.valueChanged.connect(self.update_freq_label)
        
        # Search fields
        self.gene_search.returnPressed.connect(self.apply_filters)
        self.variant_search.returnPressed.connect(self.apply_filters)
        self.region_search.returnPressed.connect(self.apply_filters)
    
    def open_vcf_file(self):
        """Open a VCF file dialog and load the selected file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open VCF File", "", "VCF Files (*.vcf *.vcf.gz)"
        )
        
        if file_path:
            self.load_vcf_file(file_path)
    
    def add_vcf_file(self):
        """Add a VCF file to the multi-file comparison"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Add VCF File", "", "VCF Files (*.vcf *.vcf.gz)"
        )
        
        if file_path:
            if self.multi_loader.add_vcf(file_path):
                self.update_comparison()
                
                # Show comparison tab if we have multiple files
                if len(self.multi_loader.vcf_loaders) > 1:
                    self.tab_widget.setTabVisible(2, True)
                    self.tab_widget.setCurrentIndex(2)
    
    def update_comparison(self):
        """Update the comparison plot with current VCFs"""
        comparisons = self.multi_loader.compare_vcfs()
        sample_names = [loader.file_name for loader in self.multi_loader.vcf_loaders]
        self.comparison_plot.update_plot(comparisons, sample_names)
    
    def load_vcf_file(self, file_path: str):
        """Load and process a VCF file"""
        self.status_bar.showMessage(f"Loading {os.path.basename(file_path)}...")
        self.progress_bar.setVisible(True)
        QApplication.processEvents()  # Update UI
        
        try:
            # Load VCF file
            success = self.vcf_loader.load_vcf(file_path)
            
            if not success:
                QMessageBox.critical(self, "Error", "Failed to load VCF file")
                return
            
            # Update model with variants
            self.variant_model.set_variants(self.vcf_loader.variants, self.vcf_loader)
            
            # Update chromosome filter
            self.chrom_filter.clear()
            self.chrom_filter.addItem("All Chromosomes", None)
            
            chrom_counts = self.variant_model.get_chromosome_counts()
            for chrom in sorted(chrom_counts.keys(), 
                              key=lambda x: int(x) if x.isdigit() else float('inf')):
                self.chrom_filter.addItem(f"Chr {chrom} ({chrom_counts[chrom]} variants)", chrom)
            
            # Start annotation in background
            self.start_annotation()
            
            # Update plots
            self.update_plots()
            
            # Update status
            variant_count = self.vcf_loader.get_variant_count()
            sample_count = len(self.vcf_loader.get_sample_names())
            self.status_bar.showMessage(
                f"Loaded {variant_count} variants from {sample_count} samples"
            )
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {str(e)}")
            self.status_bar.showMessage("Error loading file")
        finally:
            self.progress_bar.setVisible(False)
    
    def start_annotation(self):
        """Start variant annotation in a background thread"""
        self.progress_bar.setVisible(True)
        self.status_bar.showMessage("Annotating variants...")
        
        self.thread = DatabaseWorker(self.vcf_loader.variants[:100], 'all')  # Limit to 100 for demo
        self.thread.progress.connect(self.progress_bar.setValue)
        self.thread.message.connect(self.status_bar.showMessage)
        self.thread.finished.connect(self.on_annotation_finished)
        self.thread.start()
    
    def on_annotation_finished(self, annotations):
        """Handle completion of variant annotation"""
        self.variant_model.variant_annotations = annotations
        self.variant_table.update_table(self.variant_model.get_filtered_variants(), annotations)
        self.progress_bar.setVisible(False)
        self.status_bar.showMessage("Annotation complete")
    
    def update_plots(self):
        """Update all plots with current data"""
        # Chromosome distribution
        chrom_counts = self.variant_model.get_chromosome_counts()
        self.chrom_plot.update_plot(chrom_counts)
        
        # Impact distribution
        impact_counts = self.variant_model.get_impact_counts()
        self.impact_plot.update_plot(impact_counts)
        
        # Quality histogram
        quality_scores = self.variant_model.get_quality_distribution()
        self.quality_hist.update_plot(quality_scores)
        
        # Frequency histogram
        frequencies = self.variant_model.get_frequency_distribution()
        self.freq_hist.update_plot(frequencies)
    
    def apply_filters(self):
        """Apply the current filters to the variant data"""
        filters = {}
        
        # Chromosome filter
        chrom_filter = self.chrom_filter.currentData()
        if chrom_filter:
            filters['chromosome'] = chrom_filter
        
        # Impact filter
        impact_filters = []
        if self.missense_check.isChecked():
            impact_filters.append('missense')
        if self.nonsense_check.isChecked():
            impact_filters.append('nonsense')
        if self.synonymous_check.isChecked():
            impact_filters.append('synonymous')
        if self.other_check.isChecked():
            impact_filters.append('other')
        
        if impact_filters:
            filters['impact'] = impact_filters
        
        # Frequency filter
        max_freq = self.freq_slider.value() / 100.0
        if max_freq < 1.0:
            filters['max_freq'] = max_freq
        
        # Search filters
        if self.gene_search.text():
            filters['gene_search'] = self.gene_search.text()
        
        if self.variant_search.text():
            filters['variant_search'] = self.variant_search.text()
        
        if self.region_search.text():
            filters['region_search'] = self.region_search.text()
        
        # Apply filters to model
        self.variant_model.apply_filters(filters)
        
        # Update table with filtered variants
        self.variant_table.update_table(
            self.variant_model.get_filtered_variants(),
            self.variant_model.variant_annotations
        )
        
        # Update status
        filtered_count = len(self.variant_model.get_filtered_variants())
        total_count = len(self.variant_model.variants)
        self.status_bar.showMessage(
            f"Showing {filtered_count} of {total_count} variants"
        )
    
    def update_freq_label(self, value: int):
        """Update the frequency filter label"""
        self.freq_label.setText(f"Max Frequency: {value}%")
    
    def on_variant_selected(self):
        """Handle selection of a variant in the table"""
        selected_items = self.variant_table.selectedItems()
        
        if not selected_items:
            self.variant_detail.update_details(None)
            return
        
        # Get the selected row (all items in selection are from same row)
        row = selected_items[0].row()
        
        # Get the variant from the model
        filtered_variants = self.variant_model.get_filtered_variants()
        if row >= len(filtered_variants):
            return
        
        variant = filtered_variants[row]
        variant_key = f"{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0]}"
        
        # Get annotation from cache
        annotation = self.variant_model.variant_annotations.get(variant_key, {})
        
        # Update details widget
        self.variant_detail.update_details(variant, annotation)
    
    def export_variants(self, format_type: str):
        """Export filtered variants to a file"""
        if not self.vcf_loader or not self.variant_model.get_filtered_variants():
            QMessageBox.warning(self, "Export Error", "No variants to export")
            return
        
        # Get save file path
        if format_type == 'vcf':
            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save VCF File", "", "VCF Files (*.vcf)"
            )
        elif format_type == 'csv':
            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save CSV File", "", "CSV Files (*.csv)"
            )
        elif format_type == 'json':
            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save JSON File", "", "JSON Files (*.json)"
            )
        else:
            return
        
        if file_path:
            # Add extension if not provided
            if not file_path.lower().endswith(f".{format_type}"):
                file_path += f".{format_type}"
            
            # Save the file
            success = self.variant_model.save_to_file(file_path, format_type)
            
            if success:
                QMessageBox.information(self, "Export Complete", 
                                       f"Variants successfully exported to {file_path}")
            else:
                QMessageBox.critical(self, "Export Error", 
                                    "Failed to export variants")


def main():
    """Application entry point"""
    app = QApplication(sys.argv)
    
    # Set application style
    app.setStyle('Fusion')
    
    # Create and show main window
    window = MainWindow()
    window.show()
    
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
