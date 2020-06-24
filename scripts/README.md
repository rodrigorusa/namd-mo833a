Este documento descreve os _scripts_ contidos neste diretório e serve como guia para reprodução dos experimentos.

### 1. Configurando o ambiente
Nesta etapa o ambiente CLAP será instalado e será configurado as chaves de acesso ao provedor _AWS_ e _login_ de acesso as máquinas. Ao final da etapa será criado os recursos _placement group_, _security group_ e um sistema _EFS_ para armazenar os binários da compilação.
1. Instale o CLAP, siga os passos descritos em <https://clap.readthedocs.io/en/latest/introduction.html#installation-guide>
1. Copie as pastas _config/_ e _groups/_ para o diretório _~/.clap_
1. Configure o provedor de acesso
    1.  No arquivo _~/.clap/configs/providers.yaml_ veja se os valores dos campos **_access_keyfile_** e **_secret_access_keyfile_** da configuração **_namd-config-us-east-1_** correspondem com suas chaves armazenadas no diretório ~/.clap/private. Caso não correspondam, renomeie-as para ficarem com mesmo nome, pois essas chaves são utilizadas em outros _scripts_.
1. Configurando o _login_
    1. No arquivo _~/.clap/configs/logins.yaml_ edite os campos **_keypair_public_file_** e **_keypair_private_file_** da configuração **_namd-login-ubuntu_** com suas chaves armazenadas no diretório _~/.clap/private_. Neste caso não é necessário manter o mesmo nome do arquivo.
    1. Neste mesmo arquivo _~/.clap/configs/logins.yaml_, altere o campo **_keypair_name_** com seu _key_pair_ de acesso.
1. Criando _placement group_, _security group_ e sistema _EFS_
    1. Entre no diretório _namd-mo833a/script_ e execute o _script_ **_setup_aws.sh_**
    **Importante:** Este _script_ vai gerar no diretório _scripts/_ dois arquivos texto, _security_group.txt_ e _efs_ip.txt_, não remova-os, pois estes são utilizados por outros _scripts_.

### 2. Compilando o NAMD
Nesta etapa será instanciado um _cluster_ para compilar o código fonte e copiar os binários gerados no sistema _EFS_ criado na etapa anterior. Ao final o _cluster_ será destruído.
1. Entre no diretório _namd-mo833a/scripts_ execute o _script_ **_compile.sh_**.

### 3. Criando _cluster_ para simulação
Nesta etapa será instanciado um _cluster_ para execução dos experimentos.
1. Entre no diretório _namd-mo833a/scripts_ e execute o _script_ **_start_cluster.sh_** _<CFG>_
Onde: 
* _<CFG>_ é a configuração do _cluster_ a ser instanciado, podendo ser:
    * _CFG-1_: _c5.large-2x_
    * _CFG-2_: _c5.large-4x_
    * _CFG-3_: _c5.large-8x_
    * _CFG-4_: _c5.large-16x_

### 4. Executando uma simulação
Nesta etapa será executado a aplicação no _cluster_ no caso de teste especificado no cluster indicado.
1. Entre no diretório _namd-mo833a/scripts_ e execute o _script_ **_run.sh_** _<CFG>_ _<TC>_ _<stop>_ _<results_folder>_ _<args>_
Onde:
* _<CFG>_ indica a configuração do _cluster_ a ser utilizado, podendo ser:
    * _CFG-1_: _c5.large-2x_
    * _CFG-2_: _c5.large-4x_
    * _CFG-3_: _c5.large-8x_
    * _CFG-4_: _c5.large-16x_
* _<TC>_ indica o caso de teste a ser executado, podendo ser:
    * _TC-1_: _ApoA1_
    * _TC-2_: _ATPase_
    * _TC-3_: _STMV_
* _<stop>_ indica se após a execução o _cluster_ deve ser destruído, _yes_ o _cluster_ será destruído e _no_ não será destruído.
* _<results_folder>_ indica o _path_ completo onde será salvos os resultados
* _<args>_ indica argumentos a serem passados na execução. Por exemplo, caso queira interromper a execução após _N_ iterações utilize o argumento _-max-pi_ _N_

2. Caso queira executar todos os experimentos realizados neste trabalho, entre no diretório _namd-mo833a/scripts_ e execute o _script_ **_run_all.sh_** <results_folder>, onde o argumento _<results_folder>_ indica o _path_ completo onde será salvos os resultados.
**Importante:** Este _script_ pode demorar para executar por completo.

### 5. Destruindo recursos criados
Nesta etapa será destruído os recursos criados na etapa 1.5, removendo o sistema _EFS_, _placement group_ e _security group_.
1. Entre no diretório _namd-mo833a/scripts_ e execute o _script_ **_clean.sh_**

### 6. Gerando sumário e gráficos
Nesta etapa será gerado um sumário e gráficos com os resultados coletados dos experimentos.
1. Para gerar o sumário entre no diretório _namd-mo833a/scripts_ e execute _python summary.py_ -f <experimental_results>
Onde:
* _<experimental_results>_ indica o _path_ completo da pasta contendo os resultados. Por exemplo: _~/namd-833a/experimental_results_
1. Para gerar os gráficos, no mesmo diretório execute _python plots.py_ -f <experimental_results>
**Importante:** Estes _scripts_ requerem os pacotes _argparse, numpy, pandas e matplotlib_ instalados.