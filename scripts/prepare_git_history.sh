#!/bin/sh
set -eu
# Запускати з кореня проєкту на власному комп'ютері.
# Історія відображає реальні ГРУПИ логічних компонентів згенерованого AI коду.
# Вона не імітує старі дати та не приписує студенту авторство AI-коду.
if [ -d .git ]; then
  echo 'Репозиторій .git уже існує. Цей скрипт призначений для чистого проєкту.'
  exit 1
fi
if ! git config --global user.name >/dev/null 2>&1 || ! git config --global user.email >/dev/null 2>&1; then
  echo 'Спочатку налаштуйте git config user.name та user.email у власному акаунті.'
  exit 1
fi
git init
git add package.json index.html README.md sample_system.json
git commit -m '[AI] docs: requirements, variant and project skeleton'
git add src/system.js
git commit -m '[AI] feat: generate diagonally dominant systems'
git add src/solvers.js
git commit -m '[AI] feat: pivoted Gaussian, Jacobi and Seidel solvers'
git add src/charts.js src/style.css src/app.js index_offline.html
git commit -m '[AI] feat: six charts, scaling and animations'
git add tests/ scripts/generate_results.mjs
git commit -m '[AI] test: validate methods, singular cases and stopping rules'
git add docs/ scripts/prepare_git_history.sh
git commit -m '[AI] docs: UML diagrams, tests and measured results'
git branch -M main
echo 'Готово. Для GitHub: git remote add origin <url> && git push -u origin main'
