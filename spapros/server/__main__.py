"""Start the spapros server: ``python -m spapros.server --host 0.0.0.0 --port 8000 --data-dir /srv/spapros``."""

import argparse


def main() -> None:
    import uvicorn

    from spapros.server.app import create_app

    parser = argparse.ArgumentParser(description="Run the spapros web backend.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--data-dir", default=None, help="Uploads and results (default: $SPAPROS_SERVER_DATA_DIR).")
    parser.add_argument("--workers", type=int, default=None, help="Jobs run at the same time (default: 1).")
    args = parser.parse_args()
    uvicorn.run(create_app(data_dir=args.data_dir, workers=args.workers), host=args.host, port=args.port)


if __name__ == "__main__":
    main()
