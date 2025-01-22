# TMDalgen

TMDalgen to projekt przygotowany w ramach pracy inżynierskiej: "Przygotowanie oprogramowania dla projektu poszukiwania metodami globalnej optymalizacji nowych, kilkuwarstwowych, dwuwymiarowych struktur z grupy dichalkogenków".

## Spis treści

- [Wymagania](#wymagania)
- [Instalacja i uruchomienie](#instalacja-i-uruchomienie)
- [Struktura projektu](#struktura-projektu)
- [Autor](#autor)

## Wymagania
- Python 3.10.12
- SIESTA 5.2.0
- Wirtualne środowisko (zalecane)

### Wymagane biblioteki
- ase 3.23.0
- numpy 1.26.4

## Instalacja i uruchomienie

Do uruchomienia obliczeń konieczne jest wcześniejsze zainstalowanie oprogramowania SIESTA oraz dwóch zewnętrznych bibliotek. W celu uniknięcia problemów z kompatybilnością zalecane jest używanie wersji programów i bibliotek wymienionych powyżej. Zalecane jest stworzenie izolowanego, wirtualnego środowiska i instalacja bibliotek w następujący sposób: \
• python -m venv env \
• source env/bin/activate \
• pip install ase==3.23.0 \
• pip install numpy==1.26.4  

Do uruchomienia programu konieczne jest umieszczenie w katalogu projektu pliku z podstawową struktura dichalkogenka (w formacie .xyz lub .json), a także pseudopotencjałów wykorzystywanych do obliczeń dla poszczególnych pierwiastków (w formacie .psf) w folderze pseudos. Przed uruchomieniem obliczeń należy aktywować wirtualne środowisko: \
• source env/bin/activate \
oraz ustawić ścieżki do kalkulatora SIESTA i do plików pseudopotencjałów: \
• export ASE_SIESTA_COMMAND="mpirun -np liczba_rdzeni sciezka_do_folderu_siesta < PREFIX.fdf > PREFIX.out" \
• export SIESTA_PP_PATH=sciezka_do_folderu_pseudos

Przed uruchomieniem programu należy ustawić parametry algorytmu w pliku input.txt:
```
n_generations = 10              # Liczba pokoleń
pop_size = 10                   # Rozmiar populacji
n_best = 1                      # Liczba najlepszych osobników z poprzedniego pokolenia trafiająca bezpośrednio do nowego pokolenia
n_child = 2                     # Liczba nowych osobników tworzona poprzez krzyżowanie
n_mut = 2                       # Liczba nowych osobników tworzona poprzez mutację
struct_filename = MoS2.xyz      # Nazwa pliku z podstawową strukturą dichalkogenka
size = 4x4                      # Rozmiar struktur w populacji
n_atoms = 4                     # Liczba atomów pomiędzy warstwami
n_change = 2                    # Liczba atomów wymieniana pomiędzy strukturami w trakcie krzyżowania
atom_symbol = Mo                # Symbol chemiczny atomów pomiędzy warstwami
mag_moment = 1.0                # Początkowy moment magnetyczny nadawany atomom pomiędzy warstwami
label = MoS2                    # Etykieta nadawana plikom wyjściowym kalkulatora SIESTA
```
Powyżej zostały przedstawione przykładowe parametry algorytmu do obliczeń dwuwarstwowych struktur dwusiarczku molibdenu, będących powiększoną czterokrotnie w kierunku x i y komórką elementarną MoS2 z czterama atomami molibdenu umieszczonymi pomiędzy warstwami.

Po wykonaniu powyższych kroków można uruchomić progam komendą (znajdując się w katalogu projektu TMDalgen): \
• python3 main.py


W wyniku działania programu zostaną utworzone katalogi pop_ dla każdego wygenerowanego pokolenia (np. pop0 to populacja początkowa) zawierające pliki wyjściowe kalkulatora każdej analizowej struktury w oddzielnym folderze (cand_ - losowe struktury, child_ - struktury powstałe w wyniku krzyżowania, mut_ - struktury powstałe w wyniku mutacji). W głównym folderze projektu zostaną także zapisane pliki z wygnerowanymi strukturami (sorted_pop_.traj) oraz pliki zawierające energie struktur (energy_pop_.txt) dla każdego pokolenia.  

## Struktura projektu
```
TMDalgen/
├── main.py 			# Główny plik projektu
├── input.txt 			# Plik konfiguracyjny
├── MoS2.xyz 			# Plik z podstawową strukturą dichalkogenka
├── functions/ 			# Folder z modułami zawierającymi funkcje
│   ├── calculator.py 		# Funkcja tworząca kalkulator SIESTA o zadanych parametrach
│   ├── continue_generation.py 	# Funkcja do kontynuowania niezakończonego generowania pokolenia 
│   ├── crossover.py 		# Funkcja przeprowadzająca operacje krzyżowania między dwiema strukturami
│   ├── gen_energy_file.py 	# Funkcja zapisująca energie struktur w danym pokoleniu do pliku .txt
│   ├── gen_rand_struct.py 	# Funkcja generująca dwuwarstwową strukturę z losowo rozmieszczonymi atomami
│   ├── gen_rand_pop.py 	# Funkcja generująca losową populacje struktur
│   ├── load_config.py 		# Funkcja wczytująca parametry algorytmu z pliku input.txt
│   ├── mutation.py 		# Funkcja przeprowadzająca operacje mutacji struktury
│   ├── prep_generation.py 	# Funkcja przygotowująca nową populację na podstawie poprzedniego pokolenia
│   ├── prep_struct.py 		# Funkcja generująca dwuwarstwową strukturę na podstawie pliku .xyz
│   ├── small_functions.py 	# Moduł zawierający funkcje pomocnicze
│   └── sort_population.py 	# Funkcja sortująca struktury w danym pokoleniu (od najniższej do najwyższej energii)
├── pseudos/		      	# Folder z pseudopotencjałami wykorzystywanymi do obliczeń
├── docs/                 	# Dokumentacja projektu
└── README.md             	# Opis projektu
```
Szczegółowa dokumentacja projektu znajduje się w pliku docs/index.html.

## Autor
Adam Rek 2025
