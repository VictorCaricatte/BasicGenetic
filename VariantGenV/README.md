# VariantGenV: Ferramenta de Visualização e Análise de Variantes Genéticas

VariantGenV é uma aplicação desktop desenvolvida em Python usando PyQt6 e Pyqtgraph para visualização e análise de variantes genéticas a partir de arquivos VCF (Variant Call Format). Ela oferece funcionalidades para carregar, filtrar, anotar e comparar dados de variantes, bem como visualizar suas distribuições por cromossomo, tipo de impacto, qualidade e frequência alélica.

## Recursos

* **Carregamento de VCF:** Carrega arquivos VCF únicos ou múltiplos (.vcf e .vcf.gz).
* **Filtragem de Variantes:** Filtra variantes por:
    * Cromossomo.
    * Tipo de impacto (missense, nonsense, synonymous, etc.).
    * Frequência alélica máxima.
    * Termos de pesquisa por nome de gene, ID da variante (rsID) ou região genômica (por exemplo, chr1:10000-20000).
* **Anotação de Variantes (Simulado):** Inclui um sistema de anotação de variantes com cache para simular consultas a bancos de dados como dbSNP, ClinVar e gnomAD, e informações de gene. A implementação real exigiria configuração de APIs/credenciais.
* **Visualizações Interativas:**
    * **Resumo:** Gráficos de barras que mostram a distribuição de variantes por cromossomo e tipo de impacto.
    * **Qualidade/Frequência:** Histogramas para visualização das distribuições de escores de qualidade e frequências alélicas. Inclui linhas de referência (por exemplo, Q30 para qualidade, 1% para frequência).
    * **Comparação:** Quando múltiplos arquivos VCF são carregados, visualiza variantes comuns e únicas entre as amostras.
    * **Tabela de Variantes:** Uma tabela dinâmica e classificável exibindo informações detalhadas da variante, com codificação de cores para status de filtro e tipo de impacto.
* **Detalhes da Variante:** Um painel dedicado exibe informações abrangentes para a variante selecionada na tabela, incluindo dados VCF brutos, informações simuladas do gene e anotações de bancos de dados.
* **Exportação de Dados:** Exporta variantes filtradas para os formatos VCF, CSV e JSON.
* **Interface do Usuário Responsiva:** As operações de carregamento e anotação de arquivos são executadas em threads de segundo plano para evitar o congelamento da UI.
* **Ajuda Integrada:** Uma guia de ajuda fornece instruções detalhadas sobre o uso da aplicação.

## Estrutura do Código

* `DatabaseWorker` (QThread): Lida com operações intensivas de banco de dados (anotação) em um thread separado para manter a UI responsiva.
* `VariantAnnotator`: Gerencia a anotação de variantes usando caches e métodos de consulta (mocked para a maioria dos DBs, com implementação BigQuery para gnomAD se as credenciais existirem).
    * `query_dbsnp`: Simula a consulta ao dbSNP.
    * `query_clinvar`: Simula a consulta ao ClinVar.
    * `query_gnomad`: Consulta gnomAD via Google BigQuery (se configurado) ou dados mockados.
    * `get_gene_info`: Simula a obtenção de informações do gene.
* `VCFLoader`: Responsável por carregar e analisar arquivos VCF usando a biblioteca `vcf` (provavelmente `PyVCF` ou similar).
* `MultiVCFLoader`: Lida com o carregamento e a comparação de múltiplos arquivos VCF, encontrando variantes comuns e únicas.
* `VariantTableModel`: Gerencia os dados da variante, incluindo filtragem e preparação de dados para exibição e plotagem.
* `ChromosomePlotWidget`: Widget Pyqtgraph para exibir contagens de variantes por cromossomo.
* `ImpactPlotWidget`: Widget Pyqtgraph para exibir a distribuição de variantes por tipo de impacto.
* `QualityHistogramWidget`: Widget Pyqtgraph para histogramas de pontuações de qualidade.
* `FrequencyHistogramWidget`: Widget Pyqtgraph para histogramas de frequência alélica.
* `ComparisonPlotWidget`: Widget Pyqtgraph para visualizar comparações entre múltiplos arquivos VCF.
* `VariantTableWidget`: QTableWidget personalizado para exibir dados de variantes, com recursos de classificação e coloração.
* `VariantDetailWidget`: QWidget para exibir informações detalhadas sobre uma variante selecionada.
* `HelpTab`: Uma guia da interface do usuário que fornece instruções e informações sobre a aplicação.
* `MainWindow` (QMainWindow): A classe principal da aplicação que configura a interface do usuário, conecta sinais/slots e orquestra as operações.

## Como Usar

### 1. Requisitos

* Python 3.x
* `PyQt6`
* `pyqtgraph`
* `numpy`
* `pandas`
* `vobject` (a biblioteca `vcf` (provavelmente `PyVCF`) é necessária para carregar arquivos VCF).

### 2. Instalação

```bash
pip install PyQt6 pyqtgraph numpy pandas PyVCF
# Se for usar BigQuery
pip install google-cloud-bigquery
