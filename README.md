# Tradeoff entre morte e desemprego na pandemia de COVID-19

Estudo da taxa de troca no tradeoff enfrentado, no início da pandemia de covid19, entre mortes e desemprego. 

As decisões governamentais pelo isolamento físico tinham que escolher entre ter:
- alto nível de isolamento com:
  - baixa quantidade de mortes; e
  - alta quantidade de desemprego.
- baixo nível de isolamento com:
  - alta quantidade de mortes; e
  - baixa quantidade de desemprego.

A partir de um modelo epidemiológico Susceptible-Exposed-Infectious-Recovered (SEIR), em python, estimou-se a taxa de mortalidade para os diferentes níveis de isolamento, para todos os estados brasileiros, suas grandes regiões, e para o país.
Com um modelo de Vetores Auto Regressivos (VAR), em R, buscou-se estimar a relação entre a taxa de isolamento e a taxa de desemprego, para as mesmas localidades.
Em R foi feita a análise entre essas duas taxas.

***

Contém web scraping de dados de isolamento/mobilidade (repositório do Google), de vacinação, infecções, mortes etc. (repositório do DataSUS) e de desemprego (IBGE/PNAD). Faz download e tratamento dos dados básicos e os salva para utilização no modelo SEIR.
- *Na pasta 'inputs_for_seir--R/' executar o script em R '10_for_seir.R', que lerá as funções em '11_funcs_for_seir.R'. Os arquivos de dados básicos baixados serão salvos na pasta 'inputs_for_seir--R/01_data/'. Os arquivos com os dados tratados, a serem usados no script python do modelo SEIR, serão salvos na pasta 'inputs_for_seir--R/02_output_for_seir/'.*

Realiza simulação com o modelo SEIR, estimando os melhores parâmetros com otimização numérica (força bruta).
- *Os parâmetros profundos do modelo, caminhos de pastas e etc são definidos no arquivo 'param.py', na pasta 'seir--phyton/'. Nesta pasta executar o script em python 'data.py' para importar os dados tratados e salvos pelo script '10_for_seir.R' (só precisará executar uma vez). Os parâmetros 'p_d', 'r_dth' e 'r_d' do modelo SEIR são otimizados para cada estado e região do Brasil através do script 'optm_param.py', que salvará arquivo com os parâmetros ótimos (só precisará executar uma vez). Por fim, executar o script em python 'run_model.py', que carregará os objetos do módulo 'seir.py' (que carrega 'param.py'), lerá os arquivos de dados e paraâmetros e fará a simulação do modelo: por estado e/ou região, com ou sem gráficos etc, de acordo com definições de variáveis booleanas iniciais.*

Faz análises dos dados com modelos ARIMA, dessazonalização com o X13-ARIMA/SEATS. Utiliza métodos não paramétricos de splines para desenhar curvas (funções) para as séries de tempo. Cria gráficos e salva tabelas para análise.
- *Na pasta 'analisys--R/' executar o script em R '20_mob_gamma.R', que carregará funções de '21_funcs_mob_gamma.R'. O script '30_analisis.R' apresenta um resumo da análise das séries de tempo e imprime diversos gráficos.*
