# Agent Instructions for Ampel-HU-astro Repository

## Overview
This document provides guidelines for AI assistants and automated agents working with the Ampel-HU-astro repository.

## Environment Setup

### Prerequisites
- Poetry (Python dependency management)
- MongoDB (for running jobs) - **Must be running locally**
- Python 3.10+

### Starting MongoDB
Before running Ampel jobs, ensure MongoDB is running:
```bash
# Start MongoDB service (macOS with Homebrew)
brew services start mongodb-community

# Verify it's running
brew services list | grep mongodb
```

### Installation
```bash
# Install dependencies
poetry install --all-extras

# Activate virtual environment
poetry shell
```

## Configuration Management

### Building Configuration
Ampel requires a compiled configuration file to run jobs. **Do not rebuild unnecessarily** as this is computationally expensive (~1-2 minutes).

```bash
# Only rebuild if absolutely necessary
poetry run ampel config build --no-provenance -out ampel_config.yaml -stop-on-errors=0
```

**Rules:**
- ✅ Rebuild ONLY when:
  - New units are added to `ampel/conf/ampel-hu-astro/unit.yml`
  - Configuration files in `ampel/conf/` are modified
  - Dependencies change
- ❌ Do NOT rebuild when:
  - Just running existing jobs
  - Only modifying test files
  - Only modifying Python code (unless it affects config)

### Checking Existing Config
```bash
# Check if config file exists
ls -la ampel_config.yaml 2>/dev/null || echo "Config not found"

# Check if config needs rebuilding (compare timestamps)
stat ampel/conf/ampel-hu-astro/unit.yml ampel_config.yaml
```

## Running Ampel Jobs

### Basic Job Execution
```bash
# Run a job with existing config
poetry run ampel job JOB_FILE.yaml -config ampel_config.yaml
```

### Common Issues
1. **Config file missing**: Use `--no-provenance` flag or rebuild config
2. **Database conflicts**: Reset with `poetry run ampel db drop` (careful!)

## Development Workflow

### Adding New Units
1. Create the unit file in appropriate directory:
   - T0: `ampel/contrib/hu/t0/`
   - T1: `ampel/contrib/hu/t1/`
   - T2: `ampel/contrib/hu/t2/`
   - T3: `ampel/contrib/hu/t3/`
   - T4: `ampel/contrib/hu/t4/`

2. Add to unit configuration:
   ```bash
   echo "- ampel.contrib.hu.tX.NewUnitName" >> ampel/conf/ampel-hu-astro/unit.yml
   ```

3. **Rebuild config** (required for new units):
   ```bash
   poetry run ampel config build --no-provenance -out ampel_config.yaml -stop-on-errors=0
   ```

### Testing Units
```bash
# Direct instantiation test
poetry run python -c "
from ampel.contrib.hu.t4.NewUnit import NewUnit
from ampel.log.AmpelLogger import AmpelLogger
logger = AmpelLogger(console=True)
unit = NewUnit(logger=logger, param1=value1)
result = unit.do()
print(result)
"
```

## Important Notes

### Logger Usage in Units
- T4 units (and other logical units) receive a logger via dependency injection
- Use `self.logger.info()`, `self.logger.debug()`, `self.logger.warn()`, `self.logger.error()`
- ⚠️ Method is `warn()` NOT `warning()` (will cause AttributeError)
- Messages are automatically prefixed with timestamp and unit name

### Database Operations
```bash
# Drop entire database (CAUTION!)
poetry run ampel db drop

# View database stats
poetry run ampel db stat

# Export data
poetry run ampel db export OUTPUT_DIR
```

### Cleanup
```bash
# Remove temporary files
rm -f ampel_config.yaml *.log

# Clear poetry cache (if needed)
poetry cache clear --all
```

## Best Practices

1. **Check for existing solutions first**: Search the codebase before implementing new features
2. **Follow naming conventions**: Match existing patterns in the codebase
3. **Use type hints**: Ampel makes extensive use of type annotations
4. **Keep imports organized**: Follow the existing import style
5. **Test incrementally**: Test small changes before integrating complex modifications

## Troubleshooting

### Module Import Errors
```bash
# Ensure you're using poetry run
poetry run python script.py

# Or activate the virtual environment
poetry shell
python script.py
```

### Missing Configuration Values
```bash
# Check what's available in config
poetry run python -c "
from ampel.config.AmpelConfig import AmpelConfig
config = AmpelConfig.load('ampel_config.yaml')
print(list(config._config.keys()))
"
```

### Version Information
```bash
# Check Ampel version
poetry run ampel --version

# Check component versions
poetry show | grep ampel
```

## Reference Files

- Main configuration: `ampel/conf/ampel-hu-astro/`
- Example configurations: `examples/`
- Unit implementations: `ampel/contrib/hu/t{0,1,2,3,4}/`
- Tests: `tests/`

## GitHub Issues
To check open issues in the Ampel-HU-astro repository:
```bash
# Requires gh CLI authentication (https://cli.github.com/)
gh auth login
gh issue list --repo AmpelAstro/Ampel-HU-astro --state open
```

## Contact
For questions about Ampel specifically, refer to:
- Official documentation: https://ampelproject.github.io/
- GitHub repository: https://github.com/AmpelProject/Ampel-HU-astro
